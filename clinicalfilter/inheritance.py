'''
Copyright (c) 2016 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import logging


class Inheritance(object):
    """ A class for checking whether the genotypes of a trio for a variant or
    variants in a gene fit an inheritance model specific to the gene and
    chromosome that the variant/s are in.
    """
    
    def __init__(self, variants, trio, known_gene, gene, cnv_regions=None):
        """ intialise the class with the variants and trio information
        
        We have an affected child, with two parents who may or may not be
        affected. We want a model where either some alleles are transmitted to
        the child, which might contribute to mendelian inheritance of a
        disorder, or the child has a de novo mutation that might contribute to
        the child's disorder.
        
        For some genes we know whether the gene is monoallelic, or biallelic, or
        other possibilities. This informs whether genotypes for a trio could be
        causal.
        
        Args:
            variants: list of TrioGenotypes variants in a gene
            trio: a Family object for the family
            known_gene: a dictionary of inheritance types for a gene known to be
                involved with developmental disorders or None.
            gene: symbol for gene (e.g. "ARID1B")
            cnv_regions: a list of (chrom, start, end, copy_number) tuples of
                genomic regions known to be involved in CNV syndromes.
        """
        
        self.variants = variants
        self.trio = trio
        self.known_gene = known_gene
        self.gene = gene
        self.cnv_regions = cnv_regions
        
        self.father_affected = None
        self.mother_affected = None
        if self.trio.has_parents():
            self.father_affected = self.trio.father.is_affected()
            self.mother_affected = self.trio.mother.is_affected()
        
        self.chrom_inheritance = self.variants[0].get_inheritance_type()
        
        # here are the inheritance modes defined in the known gene database
        if self.known_gene is None:
            self.gene_inheritance = set(["Biallelic", "Both", "Digenic", \
                "Hemizygous", "Imprinted", "Mitochondrial", "Monoallelic", \
                "Mosaic", "Uncertain", "X-linked dominant", \
                "X-linked over-dominance"])
        else:
            self.gene_inheritance = set(self.known_gene["inh"])
        
        # if a gene has an inheritance mode of "Both", make sure we will process
        # both mono and biallelic inheritance, but remove "Both"
        if "Both" in self.gene_inheritance:
            self.gene_inheritance.add("Biallelic")
            self.gene_inheritance.add("Monoallelic")
            self.gene_inheritance.remove("Both")
    
    def get_candidate_variants(self):
        """ screen for variants that might contribute to a childs disorder
        """
        
        # ignore variants on chroms that don't match the gene inheritance
        if not self.check_inheritance_mode_matches_gene_mode():
            return []
        
        # set up the lists for candidate variants
        compound_hets = []
        candidates = []
        
        for variant in self.variants:
            self.set_trio_genotypes(variant)
            
            # check against every inheritance mode for the gene
            for inheritance in self.inheritance_modes & self.gene_inheritance:
                check = self.examine_variant(variant, inheritance)
                
                group = []
                if check == "compound_het":
                    group = compound_hets
                elif check == 'single_variant':
                    group = candidates
                
                group.append((variant, (check,), (inheritance,), (self.gene,)))
            
            logging.info("{} {}:{} {}".format(self.trio.child.get_id(),
                variant.get_chrom(), variant.get_position(), self.log_string))
        
        variants = candidates + self.check_compound_hets(compound_hets)
        
        # convert the tuples in the list of passed variants to lists, so that
        # we can merge variants with multiple inheritance modes, multiple genes,
        # multiple check types
        return [ (x[0], list(x[1]), list(x[2]), list(x[3])) for x in variants ]
    
    def check_inheritance_mode_matches_gene_mode(self):
        """ make sure that the mode of inheritance for the gene makes sense
        for the mode of inheritance for the chromosome - there's no point
        looking through hemizygous genes on an autosomal chromosome, or
        biallelic genes on an X chromosome.
        """
        
        # at this point make sure we only deal with genes that have at least
        # one correct modes of inheritance for the given chromosome type
        return len(self.inheritance_modes & self.gene_inheritance) > 0
    
    def set_trio_genotypes(self, variant):
        """ sets the genotypes for the trio as Variant objects
        
        Args:
            variant: TrioGenotypes object
        """
        
        # allow for children without parents
        if self.trio.has_parents():
            self.child = variant.child
            self.mom = variant.mother
            self.dad = variant.father
        else:
            self.child = variant.child
            self.mom = None
            self.dad = None
    
    def examine_variant(self, variant, inheritance):
        """ examines a single variant for whether or not to report it
        
        Args:
            variant: TrioGenotypes object
            inheritance: inheritance mode to check ("Monoallelic", "Biallelic")
        
        Returns:
            code for whether to add the variant to the list of flagged variants
            ("single_variant"), to check if the variant could act in concert as
            a compound het ("compound_het"), or whether to ignore the variant
            ("nothing").
        """
        
        if variant.is_cnv():
            cnv_checker = CNVInheritance(self.trio, self.known_gene, self.gene, self.cnv_regions)
            check = cnv_checker.check_single_inheritance(variant)
            self.log_string = cnv_checker.log_string
            return check
        
        if not self.trio.has_parents():
            return self.check_variant_without_parents(inheritance)
        
        if self.child.is_hom_alt():
            return self.check_homozygous(inheritance)
        elif self.child.is_het():
            return self.check_heterozygous(inheritance)
        
        self.log_string = "not hom alt nor het: " + str(self.child) + " with \
            inheritance" + self.chrom_inheritance
                           
        return "nothing"
    
    def check_if_any_variant_is_cnv(self):
        """ checks if any of the variants in a gene are CNVs
        """
        
        return any([ x.is_cnv() for x in self.variants ])
    
    def check_compound_hets(self, variants):
        """ checks for compound hets within a gene
        
        Args:
            variants: list of (TrioGenotypes, check, inh, gene) tuples
        
        Returns:
            list of variants that are compatible with being compound heterozygotes
        """
        
        if len(variants) < 2:
            return []
        
        compound = set([])
        for first in variants:
            for second in variants:
                if self.is_compound_pair(first[0], second[0]):
                    compound = compound | {first, second}
        
        return list(compound)
    
    def is_compound_pair(self, first, second):
        """ determines whether two variants form a compound pair
        
        Args:
            first: TrioGenotypes object for first variant
            second: TrioGenotypes object for second variant
        
        Returns:
            true/false for whather the pair of variants could be a compound het.
        """
        
        if first == second:
            return False
        
        # we don't include compound hets where we lack parents, and both
        # variants have missense equivalent consequences, since these are very
        # unlikely to be pathogenic.
        if not self.trio.has_parents() and \
          first.child.is_missense(self.gene) and second.child.is_missense(self.gene):
            return False
        
        # some CNVs get lumped with NA "." gene values, which mean when
        # we get two CNVs under "." gene IDs, these automatically come
        # through as compound hets, even though they might be on
        # different chroms
        if first.get_genes() == ["."]:
            return False
        
        # now we have two different variants in the same gene
        self.set_trio_genotypes(first)
        mom_1, dad_1 = self.mom, self.dad
        self.set_trio_genotypes(second)
        mom_2, dad_2 = self.mom, self.dad
        
        # if either variant is a CNV, then we check the inheritance separately
        if first.is_cnv() or second.is_cnv():
            return self.check_pair_with_cnv(first, second)
        
        # assume variants in probands without parents are compound hets
        if mom_1 is None:
            return True
        
        # compound hets on the X chromosome occur when the father has a
        # nonref genotype and is affected (or both ref, but one de novo)
        if first.get_chrom() == "X" and (dad_1.is_hom_alt() or dad_2.is_hom_alt()) and \
            not self.father_affected:
            return False
        
        if (mom_1.is_hom_ref() and mom_2.is_not_ref() \
            and dad_1.is_not_ref() and dad_2.is_hom_ref()) or \
            (mom_1.is_not_ref() and mom_2.is_hom_ref() \
            and dad_1.is_hom_ref() and dad_2.is_not_ref()):
            return True
        
        return False
    
    def check_pair_with_cnv(self, first, second):
        """ check whether a pair of variants with 1+ CNV could be a compound het
        
        Args:
            first: TrioGenotypes object for first variant
            second: TrioGenotypes object for second variant
        
        Returns:
            true/false for whather the pair of variants could be a compound het.
        """
        
        # swap the variants around so the first variant has to be a CNV. Note
        # that the other variant might still be a CNV.
        if second.is_cnv():
            first, second = second, first
        
        # if one of the variants is not a CNV, then check
        if not second.is_cnv():
            # set the parental genotypes for the SNV variant
            self.set_trio_genotypes(second)
            mom_2, dad_2 = self.mom, self.dad
            
            # get the inheritance state of the CNV variant
            inh = [first.child.format["INHERITANCE"], first.child.format["CIFER_INHERITANCE"]]
            paternal = any([ "paternal" in x for x in inh ])
            maternal = any([ "maternal" in x for x in inh ])
            
            # If the CNV is paternally inherited, then for the other variant, we
            # need it to be inherited from the mother, and not from the father.
            # This is vice-versa if the CNV is maternally inherited.
            if paternal:
                return dad_2.is_hom_ref() and mom_2.is_not_ref()
            elif maternal:
                return dad_2.is_not_ref() and mom_2.is_hom_ref()
        elif first.is_cnv() and second.is_cnv():
            return True
        
        return False

class Autosomal(Inheritance):
    
    def __init__(self, variants, trio, known_genes, gene, cnv_regions=None):
        
        super(Autosomal, self).__init__(variants, trio, known_genes, gene, cnv_regions)
        
        self.inheritance_modes = set(["Monoallelic", "Biallelic", "Both", 'Imprinted'])
    
    def check_variant_without_parents(self, inheritance):
        """ test variants in children where we lack parental genotypes
        """
        
        self.log_string = "autosomal without parents"
        if self.child.is_het() and inheritance == "Biallelic":
            return "compound_het"
        elif self.child.is_hom_alt() and inheritance == "Biallelic":
            return "single_variant"
        elif self.child.is_het() and inheritance == "Monoallelic":
            return "single_variant"
        elif self.child.is_het() and inheritance == "Imprinted" and \
                self.child.is_lof() and self.known_gene is not None and \
                'Imprinted' in self.known_gene['inh']:
            return "single_variant"
        
        return "nothing"
    
    def check_heterozygous(self, inheritance):
        """ checks if a heterozygous genotype could contribute to disease
        """
        
        if "Monoallelic" == inheritance:
            # dominant, should be reported
            report = "single_variant"
        elif "Biallelic" == inheritance:
            # recessive: should be marked for compound-het screen
            report = "compound_het"
        elif inheritance == 'Imprinted':
            report = 'nothing'
            if (self.dad.is_not_ref() or self.mom.is_not_ref()) and \
                    self.child.is_lof() and self.known_gene is not None and \
                    'Imprinted' in self.known_gene['inh']:
                self.log_string = 'possible imprinted variant'
                return "single_variant"
        
        if self.mom.is_hom_ref() and self.dad.is_hom_ref():
            self.log_string = "de novo as {0}".format(report)
            return report
        elif (self.dad.is_not_ref() and self.father_affected) and \
             (self.mom.is_hom_ref() or self.mother_affected) or \
             (self.mom.is_not_ref() and self.mother_affected) and \
             (self.dad.is_hom_ref() or self.father_affected):
            self.log_string = "transmitted from aff, other parent non-carrier or aff"
            return report
        elif "Biallelic" == inheritance and \
             ((self.dad.is_not_alt() or self.father_affected) and \
             (self.mom.is_not_alt() or self.mother_affected)):
            # TODO: could simplify this check
            self.log_string = "het-check for recessive genes and unaff parents not homoz"
            return report
        else:
            self.log_string = "typically for trios with non-de novo unaffected parents"
            return "nothing"
        
    def check_homozygous(self, inheritance):
        """ checks if a homozygous genotype could contribute to disease
        """
        
        if self.dad.is_hom_ref() or self.mom.is_hom_ref():
            # some hom alts might occur as CNV DELs change a het call to a hom,
            # catch these if a non-mendelian hom alt overlaps a CNV
            if self.check_if_any_variant_is_cnv():
                self.log_string = "non-mendelian, but CNV might affect call"
                return "compound_het"
            else:
                #NB: will miss one inherited copy and one de-novo at same site
                self.log_string = "non-mendelian trio"
                return "nothing"
        elif "Biallelic" == inheritance:
            if self.mom.is_het() and self.dad.is_het():
                self.log_string = "both parents het in biallelic gene"
                return "single_variant"
            elif (((self.mom.is_hom_alt() and self.mother_affected) and \
                (self.dad.is_not_alt() or self.father_affected)) or \
                ((self.dad.is_hom_alt() and self.father_affected) and \
                (self.mom.is_not_alt() or self.mother_affected))):
                self.log_string = "homoz parent aff"
                return "single_variant"
        elif "Monoallelic" == inheritance:
            if self.father_affected and self.mother_affected:
                self.log_string = "transmitted from affected parents"
                return "single_variant"
        elif "Imprinted" == inheritance and self.child.is_lof() and \
                self.mom.is_not_ref() and self.dad.is_not_ref() and \
                self.known_gene is not None and 'Imprinted' in self.known_gene['inh']:
            self.log_string = "imprinted variant"
            return "single_variant"
        
        self.log_string = "non-causal homozygous variant"
        return "nothing"


class Allosomal(Inheritance):
    
    def __init__(self, variants, trio, known_genes, gene, cnv_regions=None):
        
        super(Allosomal, self).__init__(variants, trio, known_genes, gene, cnv_regions)
        
        self.inheritance_modes = set(["X-linked dominant", "Hemizygous", \
            "Monoallelic", "X-linked over-dominance"])
        
        # on the X chrom, treat monoallelic and X-linked dominant modes of
        # inheritance the same
        if "Monoallelic" in self.gene_inheritance:
            self.gene_inheritance.add("X-linked dominant")
            self.gene_inheritance.remove("Monoallelic")
        
        if "X-linked over-dominance" in self.gene_inheritance:
            self.gene_inheritance.add("X-linked dominant")
            # self.gene_inheritance.remove("X-linked over-dominance")
    
    def check_variant_without_parents(self, inheritance):
        """ test variants in children where we lack parental genotypes
        """
        
        self.log_string = "allosomal without parents"
        if inheritance == "X-linked dominant":
            return "single_variant"
        elif inheritance == "Hemizygous":
            if self.child.is_het():
                return "hemizygous"
            return "single_variant"
        
        return "nothing"
    
    def check_heterozygous(self, inheritance):
        """ checks if a heterozygous genotype could contribute to disease
        """
        
        if inheritance in ["X-linked dominant", "X-linked over-dominance"]:
            report = "single_variant"
        elif "Hemizygous" == inheritance:
            # recessive: should be marked for compound-het screen
            report = "compound_het"
        else:
            raise ValueError("unknown gene inheritance: " + str(inheritance))
        
        if self.mom.is_hom_ref() and self.dad.is_hom_ref():
            self.log_string = "female x chrom de novo"
            return "single_variant"
        elif (self.dad.is_hom_alt() and self.father_affected) and \
             (self.mom.is_hom_ref() or self.mother_affected) or \
             (self.mom.is_not_ref() and self.mother_affected) and \
             (self.dad.is_hom_ref() or self.father_affected):
            self.log_string = "x chrom transmitted from aff, other parent non-carrier or aff"
            return report
        elif inheritance == "X-linked over-dominance" and self.dad.is_not_ref() \
                and self.child.is_lof() and self.known_gene is not None and \
                'X-linked over-dominance' in self.known_gene['inh']:
            self.log_string = "variant inherited from dad in X-linked over-dominance gene"
            return report
        else:
            self.log_string = "variant not compatible with being causal"
            return "nothing"
        
    def check_homozygous(self, inheritance):
        """ checks if a homozygous genotype could contribute to disease
        """
        
        if inheritance not in ["X-linked dominant", "Hemizygous", "X-linked over-dominance"]:
            raise ValueError("unknown gene inheritance: " + str(inheritance))
        
        # treat male sex inheritance differently from female sex inheritance
        if self.trio.child.is_male():
            if self.mom.is_hom_ref():
                self.log_string = "male X chrom de novo"
                return "single_variant"
            elif inheritance == "X-linked over-dominance" and self.mom.is_het() \
                    and not self.mother_affected:
                self.log_string = "inherited variant in X-linked over-dominant gene, but from unaffected mother"
                return "nothing"
            elif (self.mom.is_het() and not self.mother_affected) or \
                 (self.mom.is_hom_alt() and self.mother_affected):
                self.log_string = "male X chrom inherited from het mother or hom affected mother"
                return "single_variant"
        else:
            if self.dad.is_hom_ref() or self.mom.is_hom_ref():
                # some hom alts might occur as CNV DELs change a het call to a hom,
                # catch these if a non-mendelian hom alt overlaps a CNV
                if self.check_if_any_variant_is_cnv():
                    self.log_string = "non-mendelian, but CNV might affect call"
                    return "compound_het"
                else:
                    self.log_string = "non-mendelian trio"
                    return "nothing"
            elif (self.mom.is_het() or \
                 (self.mom.is_hom_alt() and self.mother_affected)) and \
                 (self.dad.is_hom_alt() and self.father_affected):
                self.log_string = "testing"
                return "single_variant"
        
        self.log_string = "variant not compatible with being causal"
        return "nothing"
   

class CNVInheritance(object):
    
    def __init__(self, trio, known_gene, gene, cnv_regions):
        """ intialise the class
        
        Args:
            trio: family trio object
            known_genes: dictionary of known genes, currently the DDG2P set
        """
        
        self.trio = trio
        self.known_gene = known_gene
        self.gene = gene
        self.cnv_regions = cnv_regions
    
    def check_single_inheritance(self, variant):
        """ checks if a CNV could be causal by itself
        
        Args:
            variant: TrioGenotypes object for the CNV.
        
        Returns:
            "single_variant", "compound_het", or "nothing" depending on whether
            the variant could possibly contribute to the childs disorder.
        """
        
        if not self.trio.has_parents():
            return self.check_variant_without_parents(variant)
        
        # check that the inheritance status is consistent with the parental
        # affected status
        inh = [variant.child.format["INHERITANCE"], variant.child.format["CIFER_INHERITANCE"]]
        if not self.inheritance_matches_parental_affected_status(variant, inh):
            if self.check_compound_inheritance(variant):
                self.log_string = "possible compound het CNV"
                return "compound_het"
            self.log_string = "not consistent with parental affected status"
            return "nothing"
        
        if self.passes_nonddg2p_filter(variant):
            return "single_variant"
        elif self.known_gene is not None and self.passes_ddg2p_filter(variant):
            return "single_variant"
        elif self.cnv_regions is not None and self.check_cnv_region_overlap(variant, self.cnv_regions):
            return "single_variant"
        
        if self.check_compound_inheritance(variant):
            self.log_string = "possible compound het CNV"
            return "compound_het"
        
        return "nothing"
    
    def check_compound_inheritance(self, variant):
        """ checks if a CNV could contribute to a compound het
        
        Compound CNVs can occur with CNVs with copy number of 1 or 3, and if in
        a DDG2P gene, the disorder relating to the gene must be inherited in a
        biallelic or hemizygous mode. For non-compound filtering, CNVs are
        checked to see if their inheritance state (paternal, maternal) matches
        the parents affected status (maternally affected requires. In contrast,
        compound CNVs do not require inheritance = affected status, as the
        compound variant is typically incomplete in the parents, and therefore
        is not expected to alter their affected status.
        
        Candidates for compound CNVs undergo similar filtering to single CNVs.
          - CNVs covering a DDG2P gene are included if the DDG2P gene is
            biallelic or hemizygous. In comparison to the single variant
            filtering, compound biallelic CNVs can have copy number 1 or 3,
            rather than requiring a copy number of 0. If the DDG2P gene is
            hemizygous, the CNV has to be on chr X, in a female proband, and
            with a copy number status of 1 (since CN = 3 is captured as single
            variant).
          - CNVs not in DDG2P genes are included if they span > 500000 bp.
        
        Candidate compound CNVs are checked against all of the genes that they
        span, so that if they overlap another candidate compound variant (CNV
        or SNV), the variants are included in the clinical filtering output.
        
        Note that CNVs might change the aparent state of SNVs within their
        boundaries, so that a deletion might alter a heterozygous SNV to a
        homozygous alternate allele genotype. We make some checks of hom alt
        SNVs to see if their gene includes a CNV variant.
        
        Args:
            variant: TrioGenotypes object for the CNV.
        """
        
        # we don't want CNVs that don't have copy number of 1 or 3, since
        # copy number = 1 or 3 are the only ones that could operate as compound
        # hets (other copy numbers such as 0 are implicitly dominant)
        if variant.child.info["CNS"] not in {"1", "3"}:
            return False
        
        # for compound hets, we don't have to worry about checking whether the
        # inheritance ststaus is consistent with the parental affected status
        # and in fact, most compound hets would not have affected parents
        if self.passes_nonddg2p_filter(variant):
            return True
        
        # now check the DDG2P genes. We only want CNVs in genes with Biallelic
        # (copy number = 1 or 3) or Hemizygous (copy number = 1) inheritance.
        # TODO: I might be including some CNVs erroneously here, since if a SNV
        # is on a biallelic gene, but the CNV spans multiple genes, if any of
        # those genes involves a disorder inherited in a biallelic mode, then
        # the CNV will be passed through for compound het checking.
        if self.known_gene is not None:
            if "Biallelic" in self.known_gene["inh"] or \
                 (variant.get_chrom() == "X" and \
                 "Hemizygous" in self.known_gene["inh"] and \
                 not self.trio.child.is_male() and variant.child.info["CNS"] == "1"):
                return True
        
        return False
    
    def check_variant_without_parents(self, variant):
        """ check for CNVs, without relying upon any parental genotypes
        
        Args:
            variant: TrioGenotypes object for the CNV.
        """
        
        # make sure the variant has an inheritance state of "unknown" for
        # the passes_non_ddg2p_filter()
        variant.child.format["INHERITANCE"] = "unknown"
        
        if self.passes_nonddg2p_filter(variant):
            return "single_variant"
        elif self.known_gene is not None and self.passes_ddg2p_filter(variant):
            return "single_variant"
        
        if self.check_compound_inheritance(variant):
            self.log_string = "possible compound het CNV"
            return "compound_het"
        
        return "nothing"
    
    def inheritance_matches_parental_affected_status(self, variant, inh):
        """ check that the inheritance matches the parental affected status.
        
        If the variant has been inherited from the mother (ie maternally), we
        expect the mother to also be affected. For some variants we don't know
        whether how the variant was transmitted (due to uncertainties in
        classifying the transmission). In this case, we assume the inheritance
        state would be correct for the parental affected states.
        
        Args:
            variant: TrioGenotypes object for the CNV.
            inh: list of inheritance statuses of a CNV. We have two inheritance
                classifications, from VICAR (classified from parental likelihoods
                from array CGH data) and CIFER (classified from exome based read
                depths in populations). This gives lists such as:
                [maternal, maternal_inh], [not_inherited, deNovo] etc
            
        Returns:
            True/False for whether the inheritance is consistent with the
               parental affected statuses
        """
        
        # figure out whether the inheritance classifications indicate whether
        # the variant is paternally, maternally, or biparentally inherited
        paternal = any(["paternal" in x for x in inh])
        maternal = any(["maternal" in x for x in inh])
        biparental = any([y in x for x in inh for y in ["biparental", "inheritedDuo"]])
        
        if not (paternal or maternal or biparental):
            # if the variant isn't inherited (or the inheritance isn't known),
            # then the parental affected statuses are irrelevant.
            return True
        elif (paternal and self.trio.father.is_affected()) or \
              (maternal and self.trio.mother.is_affected()) or \
              (biparental and int(variant.child.info["CNS"]) == 0) or \
              (biparental and \
              (self.trio.father.is_affected() or self.trio.mother.is_affected())) or \
              (self.trio.child.is_male() and maternal and variant.get_chrom() == "X" and not self.trio.mother.is_affected()) :
            # if the inheritance status indiates that the CNV was inherited,
            # the pertinent parents need to be also affected.
            return True
        
        return False
    
    def passes_nonddg2p_filter(self, variant):
        """ checks if a CNV passes the non DDG2P criteria
        
        Args:
            variant: TrioGenotypes object for the CNV.
        """
        
        inh = variant.child.format["INHERITANCE"]
        geno = variant.child.genotype
        
        # CNVs not in known genes are check for their length. Longer CNVs are
        # more likely to be disruptively causal, and non-artifacts. The length
        # required depends on whether the CNV was inherited, and whether the
        # CNV is a deletion, or duplication
        min_len = 1000000
        
        # reportable CNVs must be longer than the minimum length
        if float(variant.child.info["SVLEN"]) >= min_len:
            self.log_string = "non-DDG2P " + geno + " CNV, inh:" + inh
            return True
        
        self.log_string = "short non-DDG2P " + geno + " CNV, inh:" + inh
        return False
    
    def passes_ddg2p_filter(self, variant):
        """ checks if a CNV passes the DDG2P CNV criteria
        
        Args:
            variant: TrioGenotypes object for the CNV.
        """
        
        if self.known_gene is None:
            self.log_string = "non-reported CNV"
            return False
        
        self.log_string = "non-reported DDG2P CNV"
        gene_type = self.known_gene["status"]
        
        if "both dd and if" in gene_type:
            self.log_string = "Both DD and IF DDG2P gene"
            return True
        elif {"confirmed dd gene", "probable dd gene"} & gene_type == set():
            return False
        
        for inh in self.known_gene["inh"]:
            if self.passes_gene_inheritance(variant, inh) or \
                    self.passes_intragenic_dup(variant, inh):
                self.log_string = "DDG2P CNV"
                return True
        
        return False
    
    def passes_gene_inheritance(self, variant, inh):
        """ create CNV filter values dependent on DDG2P inheritance mode
        
        Args:
            variant: TrioGenotypes object for the CNV.
            gene: gene ID (eg ATRX), which is in the self.known_gene dictionary
            inh: inheritance state for the gene (eg Biallelic, Monoallelic)
        
        Returns:
            True/False for whether the gene and inheritance are consistent with
            the copy number and mechanism.
        """
        
        copies = {"0", "1", "3"}
        chroms = set([str(x) for x in range(1, 23)]) | {"X"}
        mechanisms = {"Uncertain", "Loss of function", "Dominant negative", \
            "Increased gene dosage"}
        
        if inh == "Biallelic":
            copies = {"0"}
            mechanisms = {"Uncertain", "Loss of function", "Dominant negative"}
        elif inh == "Monoallelic":
            pass
        elif inh == "X-linked dominant":
            chroms = {"X"}
        elif inh == "Hemizygous":
            chroms = {"X"}
            if not variant.child.is_male():
                copies = {"3"}
                mechanisms = {"Increased gene dosage"}
        else:
            # exclude other inheritance modes (Mosaic etc) with impossible
            # criteria
            copies = {"XXXX"}
        
        # check if the CNV is a duplication that surrounds a gene, where the
        # mechanism is loss of function, since these whole-gene duplications
        # won't disrupt the gene.
        start, end = variant.child.get_range()
        surrounding_disruptive_dup = variant.child.genotype == "DUP" and \
            "Loss of function" in self.known_gene["inh"][inh] and \
            inh in ["Monoallelic", "Hemizygous", "X-linked dominant"] and \
            start < self.known_gene["start"] and end > self.known_gene["end"]
        
        return variant.get_chrom() in chroms and \
            variant.child.info["CNS"] in copies and \
            len(self.known_gene["inh"][inh] & mechanisms) > 0 and \
            not surrounding_disruptive_dup
    
    def passes_intragenic_dup(self, variant, inh):
        """ checks if the CNV is an intragenic dup (in an appropriate gene)
        
        Args:
            variant: TrioGenotypes object for the CNV.
            gene: symbol for the gene (e.g. "ARID1B")
            inh: inhritance mode for the gene (e.g. "Monoallelelic")
        """
        
        allowed_inh = set(["Monoallelic", "X-linked dominant"])
        allowed_mech = set(["Loss of Function"])
        
        # only allow duplications in dominant genes
        if variant.child.genotype != "DUP" or inh not in allowed_inh:
            return False
        
        # only allow genes with loss-of-function mechanisms
        if len(self.known_gene["inh"][inh] & allowed_mech) < 1:
            return False
        
        # find the CNV start and end, as well as the gene start and end
        cnv_start, cnv_end = variant.get_range()
        
        gene_start = self.known_gene["start"]
        gene_end = self.known_gene["end"]
        
        # check if any part of the gene is outside the CNV boundaries
        return gene_start < cnv_start or gene_end > cnv_end
    
    def check_cnv_region_overlap(self, variant, cnv_regions):
        """ finds CNVs that overlap DECIPHER syndrome regions
        
        We have a set of genome regions known to be involved in CNV disorders.
        We want to check if the current CNV has enough overlap with any of those
        regions. This function checks the positions of the CNV and the genome
        region, to determine how much the CNV overlaps the genome-region, and
        how much the genome-region overlaps the CNV. We want the CNV to have
        high overlap of the genome-region in order to claim that the CNV might
        contribute to the probands disorder.
        
        Args:
            variant: TrioGenotypes object for the CNV.
            cnv_regions: a list of (chrom, start, end, copy_number) tuples of
                genomic regions known to be involved in CNV syndromes.
        
        Returns:
            true/false for whether the current CNV overlaps any of the syndrome
            regions.
        """
        
        chrom = variant.child.get_chrom()
        start, end = variant.child.get_range()
        copy_number = int(variant.child.info["CNS"])
        
        for region_key in cnv_regions:
            region_chrom = region_key[0]
            region_start = int(region_key[1])
            region_end = int(region_key[2])
            region_copy_number = int(cnv_regions[region_key])
            
            if region_chrom != chrom or copy_number != region_copy_number:
                continue
            
            if start <= region_end and end >= region_start and \
                    self.has_enough_overlap(start, end, region_start, region_end):
                self.log_string = "in DECIPHER syndrome region"
                return True
        
        return False
    
    def has_enough_overlap(self, start, end, region_start, region_end):
        """ finds if a CNV and another chrom region share sufficient overlap
        
        Args:
            start: start position of the CNV
            end: end position of the CNV
            region_start: start position of the genome region
            region_end: end position of the genome region
        
        Returns:
            true/false for whether the CNV has sufficient overlap of the genome
            region.
        """
        
        # find the point where the overlap starts
        overlap_start = region_start
        if region_start <= start <= region_end:
            overlap_start = start
        
        # find the point where the overlap ends
        overlap_end = region_end
        if region_start <= end <= region_end:
            overlap_end = end
        
        distance = (overlap_end - overlap_start) + 1
        
        # adjust the positions before we try to divide by zero, if the "region"
        # is actually a SNV
        if end == start:
            start -= 1
        if region_end == region_start:
            region_start -= 1
        
        forward = float(distance)/(abs(end - start) + 1)
        reverse = float(distance)/(abs(region_end - region_start) + 1)
        
        # determine whether there is sufficient overlap
        return forward > 0 and reverse > 0.5
        
