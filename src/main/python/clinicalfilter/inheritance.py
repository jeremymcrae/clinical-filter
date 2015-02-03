""" A class for checking whether the genotypes of a trio for a variant or
variants in a gene fit an inheritance model specific to the gene and
chromosome that the variant/s are in.
"""

import logging


class Inheritance(object):
    """figure out whether trio genotypes fit mendelian inheritance models
    """
    
    def __init__(self, variants, trio, known_genes, gene_inheritance=None, cnv_regions=None):
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
            gene_inheritance: inheritance based on a database of genes known to
                be involved in disorders. If no gene list is provided, default
                to including every possible inheritance type.
        """
        
        self.variants = variants
        self.trio = trio
        self.known_genes = known_genes
        self.cnv_regions = cnv_regions
        
        if self.trio.has_parents():
            self.father_affected = self.trio.father.is_affected()
            self.mother_affected = self.trio.mother.is_affected()
        else:
            self.father_affected = None
            self.mother_affected = None
        
        self.chrom_inheritance = self.variants[0].get_inheritance_type()
        
        # here are the inheritance modes defined in the known gene database
        if gene_inheritance is None:
            self.gene_inheritance = set(["Biallelic", "Both", "Digenic", \
                "Hemizygous", "Imprinted", "Mitochondrial", "Monoallelic", \
                "Mosaic", "Uncertain", "X-linked dominant", \
                "X-linked over-dominance"])
        else:
            self.gene_inheritance = set(gene_inheritance)
        
        # if a gene has an inheritance mode of "Both", make sure we will process
        # both mono and biallelic inheritance, but remove "Both"
        if "Both" in self.gene_inheritance:
            self.gene_inheritance.add("Biallelic")
            self.gene_inheritance.add("Monoallelic")
            self.gene_inheritance.remove("Both")
        
        # set up the lists for candidate variants
        self.compound_hets = []
        self.candidates = []
    
    def get_candidate_variants(self):
        """ screen for variants that might contribute to a childs disorder
        """
        
        # ignore variants on chroms that don't match the gene inheritance
        if not self.check_inheritance_mode_matches_gene_mode():
            return []
        
        for variant in self.variants:
            self.set_trio_genotypes(variant)
            self.is_lof = variant.child.is_lof()
            
            # check against every inheritance mode for the gene
            for inheritance in self.inheritance_modes & self.gene_inheritance:
                check = self.examine_variant(variant, inheritance)
                self.add_variant_to_appropriate_list(variant, check, inheritance)
            
            logging.debug("{0} position {1} {2}".format(self.trio.child.get_id(),
                variant.get_position(), self.log_string))
        
        self.compound_hets = self.check_compound_hets(self.compound_hets)
        
        possibilities = self.candidates + self.compound_hets
        
        return possibilities
    
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
    
    def add_variant_to_appropriate_list(self, variant, check, inheritance):
        """ add processed variants to the appropriate list
        
        Args:
            variant: TrioGenotypes object (CNV or SNV)
            check: single_variant or compound_het
            inheritance: inheritance mode (eg Monoallelelic, Biallelic etc)
        """
        
        if check == "compound_het":
            self.compound_hets.append((variant, check, inheritance))
        elif check == "single_variant":
            self.candidates.append((variant, check, inheritance))
    
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
            self.log_string = "skipping CNVs for now"
            return "nothing"
            # cnv_checker = CNVInheritance(variant, self.trio, self.known_genes, self.cnv_regions)
            # check = cnv_checker.check_single_inheritance()
            # self.log_string = cnv_checker.log_string
            # return check
        
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
        
        for var in self.variants:
            if var.is_cnv():
                return True
        
        return False
    
    def check_compound_hets(self, variants):
        """ checks for compound hets within a gene
        """
        
        if len(variants) < 2:
            return []
        
        # check for proband without parental genotypes
        if not self.trio.has_parents():
            return variants
        
        compound = set([])
        for first in variants:
            for second in variants:
                if first[0] == second[0]:
                    continue
                
                # some CNVs get lumped with NA "." gene values, which mean when
                # we get two CNVs under "." gene IDs, these automatically come
                # through as compound hets, even though they might be on
                # different chroms
                if first[0].get_gene() == ".":
                    continue
                
                # include CNVs in compound hets
                if first[0].is_cnv() or second[0].is_cnv():
                    compound = compound | {first, second}
                    continue
                
                # now we have two different variants in the same gene
                self.set_trio_genotypes(first[0])
                mom_1 = self.mom
                dad_1 = self.dad
                self.set_trio_genotypes(second[0])
                mom_2 = self.mom
                dad_2 = self.dad
                
                # compound hets on the X chromosome occur when the father has a
                # nonref genotype and is affected (or both ref, but one de novo)
                if first[0].child.get_inheritance_type() != "autosomal" and \
                   (dad_1.is_hom_alt() or dad_2.is_hom_alt()):
                    if not self.father_affected:
                        continue
                
                # check for 111, 111 combo
                if mom_1.is_not_ref() and mom_2.is_not_ref() and \
                    dad_1.is_not_ref() and dad_2.is_not_ref():
                    # if both variants are 1/1/1, both parents must be affected
                    if self.mother_affected and self.father_affected:
                        compound = compound | {first, second}
                elif (mom_1.is_hom_ref() and dad_1.is_hom_ref()) or \
                    (mom_2.is_hom_ref() and dad_2.is_hom_ref()):
                    # one is de novo, so they both definitely get reported
                    compound = compound | {first, second}
                elif not ((mom_1.is_hom_ref() and mom_2.is_hom_ref()) or \
                    (dad_1.is_hom_ref() and dad_2.is_hom_ref())):
                    compound = compound | {first, second}
        
        return list(compound)


class Autosomal(Inheritance):
    
    def __init__(self, variants, trio, known_genes, gene_inheritance=None, cnv_regions=None):
        
        super(Autosomal, self).__init__(variants, trio, known_genes, gene_inheritance, cnv_regions)
        
        self.inheritance_modes = set(["Monoallelic", "Biallelic", "Both"])
    
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
        
        if self.mom.is_hom_ref() and self.dad.is_hom_ref():
            self.log_string = "de novo as {0}".format(report)
            return report
        # elif self.is_lof and not self.father_affected and not \
        #         self.mother_affected and inheritance == "Monoallelic" and \
        #         inheritance in self.known_genes[self.child.gene] and \
        #         self.known_genes[self.child.gene][inheritance] != "Loss of Function" and \
        #         (self.dad.is_het() or self.mom.is_het()):
        #     self.log_string = "reduced penetrance transmitted from unaffected parents"
        #     return report
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
            # elif self.is_lof and not self.father_affected and not \
            #         self.mother_affected and self.dad.is_not_ref() and \
            #         inheritance == "Monoallelic" and \
            #         inheritance in self.known_genes[self.child.gene] and \
            #         self.known_genes[self.child.gene][inheritance] != "Loss of Function" and \
            #         self.mom.is_not_ref():
            #     self.log_string = "reduced penetrance transmitted from unaffected parents"
            #     return "single_variant"
        
        self.log_string = "non-causal homozygous variant"
        return "nothing"


class Allosomal(Inheritance):
    
    def __init__(self, variants, trio, known_genes, gene_inheritance=None, cnv_regions=None):
        
        super(Allosomal, self).__init__(variants, trio, known_genes, gene_inheritance, cnv_regions)
        
        # self.inheritance_modes = set(["X-linked dominant", "Hemizygous", \
        #     "Monoallelic", "X-linked over-dominance"])
        self.inheritance_modes = set(["X-linked dominant", "Hemizygous", \
            "Monoallelic"])
        
        # on the X chrom, treat monoallelic and X-linked dominant modes of
        # inheritance the same
        if "Monoallelic" in self.gene_inheritance:
            self.gene_inheritance.add("X-linked dominant")
            self.gene_inheritance.remove("Monoallelic")
    
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
        
        if "X-linked dominant" == inheritance:
            report = "single_variant"
        elif "Hemizygous" == inheritance:
            # recessive: should be marked for compound-het screen
            report = "compound_het"
        else:
            raise ValueError("unknown gene inheritance: " + str(inheritance))
        
        if self.mom.is_hom_ref() and self.dad.is_hom_ref():
            self.log_string = "female x chrom de novo"
            return "single_variant"
        # elif self.is_lof and not self.father_affected and not \
        #         self.mother_affected and inheritance == "X-linked dominant" and \
        #         inheritance in self.known_genes[self.child.gene] and \
        #         self.known_genes[self.child.gene][inheritance] != "Loss of Function" and \
        #         self.dad.is_hom_ref() and self.mom.is_het():
        #     self.log_string = "reduced penetrance transmitted from unaffected mother"
        #     return report
        elif (self.dad.is_hom_alt() and self.father_affected) and \
             (self.mom.is_hom_ref() or self.mother_affected) or \
             (self.mom.is_not_ref() and self.mother_affected) and \
             (self.dad.is_hom_ref() or self.father_affected):
            self.log_string = "x chrom transmitted from aff, other parent non-carrier or aff"
            return report
        # elif inheritance == "X-linked over-dominance" and not \
        #       self.father_affected and \
        #      ((self.mom.is_hom_ref() and self.dad.is_hom_alt()) or \
        #      self.mom.is_het() and self.mother_affected):
        #     self.log_string = "X-linked inheritance with unaffected hom alt \
        #         males and females, but affected het females (eg PCDH19)"
        #     return report
        else:
            self.log_string = "variant not compatible with being causal"
            return "nothing"
        
    def check_homozygous(self, inheritance):
        """ checks if a homozygous genotype could contribute to disease
        """
        
        if inheritance not in ["X-linked dominant", "Hemizygous"]:
            raise ValueError("unknown gene inheritance: " + str(inheritance))
        
        # treat male sex inheritance differently from female sex inheritance
        if self.trio.child.is_male():
            if self.mom.is_hom_ref():
                self.log_string = "male X chrom de novo"
                return "single_variant"
            elif (self.mom.is_het() and not self.mother_affected) or \
                 (self.mom.is_hom_alt() and self.mother_affected):
                self.log_string = "male X chrom inherited from het mother or hom affected mother"
                return "single_variant"
            # elif self.is_lof and not self.mother_affected and \
            #         inheritance == "X-linked dominant" and \
            #         inheritance in self.known_genes[self.child.gene] and \
            #         self.known_genes[self.child.gene][inheritance] != "Loss of Function" and \
            #         self.mom.is_not_ref():
            #     self.log_string = "reduced penetrance transmitted from unaffected mother"
            #     return "single_variant"
        
        elif self.trio.child.is_female():
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
            # elif self.is_lof and not self.father_affected and not \
            #         self.mother_affected and inheritance == "X-linked dominant" and \
            #         inheritance in self.known_genes[self.child.gene] and \
            #         self.known_genes[self.child.gene][inheritance] != "Loss of Function" and \
            #         self.dad.is_hom_ref() and self.mom.is_not_ref():
            #     self.log_string = "reduced penetrance transmitted from unaffected mother"
            #     return "single_variant"
        
        self.log_string = "variant not compatible with being causal"
        return "nothing"
   

class CNVInheritance(object):
    
    def __init__(self, variant, trio, known_genes, cnv_regions):
        """ intialise the class
        
        Args:
            variant: CNV to check for inheritance
            trio: family trio object
            known_genes: dictionary of known genes, currently the DDG2P set
        """
        
        self.variant = variant
        self.trio = trio
        self.known_genes = known_genes
        self.cnv_regions = cnv_regions
    
    def check_single_inheritance(self):
        """ checks if a CNV could be causal by itself
        
        Returns:
            "single_variant", "compound_het", or "nothing" depending on whether
            the variant could possibly contribute to the childs disorder.
        """
        
        if not self.trio.has_parents():
            self.check_variant_without_parents()
        
        # check that the inheritance status is consistent with the parental
        # affected status
        inh = self.variant.child.format["INHERITANCE"]
        if not self.inheritance_matches_parental_affected_status(inh):
            if self.check_compound_inheritance():
                self.log_string = "possible compound het CNV"
                return "compound_het"
            self.log_string = "not consistent with parental affected status"
            return "nothing"
        
        if self.passes_nonddg2p_filter():
            return "single_variant"
        elif self.known_genes is not None and self.passes_ddg2p_filter():
            return "single_variant"
        elif self.cnv_regions is not None and self.check_cnv_region_overlap(self.cnv_regions):
            return "single_variant"
        # elif self.check_compound_inheritance():
        #     return "compound_het"
        
        if self.check_compound_inheritance():
            self.log_string = "possible compound het CNV"
            return "compound_het"
        
        return "nothing"
    
    def check_compound_inheritance(self):
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
        """
        
        # return True
        
        # we don't want CNVs that don't have copy number of 1 or 3, since
        # copy number = 1 or 3 are the only ones that could operate as compound
        # hets (other copy numbers such as 0 are implicitly dominant)
        if self.variant.child.info["CNS"] not in {"1", "3"}:
            return False
        
        # for compound hets, we don't have to worry about checking whether the
        # inheritance ststaus is consistent with the parental affected status
        # and in fact, most compound hets would not have affected parents
        if self.passes_nonddg2p_filter():
            return True
        
        # now check the DDG2P genes. We only want CNVs in genes with Biallelic
        # (copy number = 1 or 3) or Hemizygous (copy number = 1) inheritance.
        # TODO: I might be including some CNVs erroneously here, since if a SNV
        # is on a biallelic gene, but the CNV spans multiple genes, if any of
        # those genes involves a disorder inherited in a biallelic mode, then
        # the CNV will be passed through for compound het checking.
        genes = self.variant.child.get_genes()
        for gene in genes:
            if self.known_genes is not None and gene in self.known_genes:
                if "Biallelic" in self.known_genes[gene]["inh"] or \
                     (self.variant.get_chrom() == "X" and \
                     "Hemizygous" in self.known_genes[gene]["inh"] and \
                     self.trio.child.is_female() and \
                     self.variant.child.info["CNS"] == "1"):
                    return True
        
        return False
    
    def check_variant_without_parents(self):
        """ check for CNVs, without relying upon any parental genotypes
        """
        
        # make sure the variant has an inheritance state of "unknown" for
        # the passes_non_ddg2p_filter()
        self.variant.child.format["INHERITANCE"] = "unknown"
        
        if self.passes_nonddg2p_filter():
            return "single_variant"
        elif self.known_genes is not None and self.passes_ddg2p_filter():
            return "single_variant"
        # elif self.check_compound_inheritance():
        #     return "compound_het"
        
        if self.check_compound_inheritance():
            self.log_string = "possible compound het CNV"
            return "compound_het"
        
        return "nothing"
    
    def inheritance_matches_parental_affected_status(self, inh):
        """ check that the inheritance matches the parental affected status
        
        Args:
            inh: inheritace status of a CNV, eg maternal, deNovo etc
            
        Returns:
            True/False for whether the inheritance is consistent with the
               parental affected statuses
        """
        
        # if the inheritance status indiates that the CNV was inherited, check
        # that the pertinent parents actually are affected.
        if inh not in ["paternal", "maternal", "biparental", "inheritedDuo"]:
            return True
        
        elif (inh == "paternal" and self.trio.father.is_affected()) or \
            (inh == "maternal" and self.trio.mother.is_affected()) or \
            ((inh == "biparental" or inh == "inheritedDuo") and \
            (self.trio.father.is_affected() or self.trio.mother.is_affected())):
            return True
        
        return False
    
    def passes_nonddg2p_filter(self):
        """ checks if a CNV passes the non DDG2P criteria
        """
        
        inh = self.variant.child.format["INHERITANCE"]
        geno = self.variant.child.genotype
        
        # CNVs not in known genes are check for their length. Longer CNVs are
        # more likely to be disruptively causal, and non-artifacts. The length
        # required depends on whether the CNV was inherited, and whether the
        # CNV is a deletion, or duplication
        min_len = 1000000
        
        # reportable CNVs must be longer than the minimum length
        if float(self.variant.child.info["SVLEN"]) >= min_len:
            self.log_string = "non-DDG2P " + geno + " CNV, inh:" + inh
            return True
        
        self.log_string = "short non-DDG2P " + geno + " CNV, inh:" + inh
        return False
    
    def passes_ddg2p_filter(self):
        """ checks if a CNV passes the DDG2P CNV criteria
        """
            
        genes = self.variant.child.get_genes()
        
        self.log_string = "non-reported CNV"
        for gene in genes:
            if self.known_genes is None or gene not in self.known_genes:
                continue
            
            self.log_string = "non-reported DDG2P CNV"
            gene_type = self.known_genes[gene]["status"]
            
            inh_passes = []
            if "Both DD and IF" in gene_type:
                self.log_string = "Both DD and IF DDG2P gene"
                inh_passes.append(True)
            elif {"Confirmed DD Gene", "Probable DD gene"} & gene_type == set():
                # drop out genes with insufficient evidence
                continue
            
            for inh in self.known_genes[gene]["inh"]:
                passes = self.passes_gene_inheritance(gene, inh) or \
                    self.passes_intragenic_dup(gene, inh)
                inh_passes.append(passes)
            
            if any(inh_passes):
                self.log_string = "DDG2P CNV"
                return True
        
        return False
    
    def passes_gene_inheritance(self, gene, inh):
        """ create CNV filter values dependent on DDG2P inheritance mode
        
        Args:
            gene: gene ID (eg ATRX), which is in the self.known_genes dictionary
            inh: inheritance state for the gene (eg Biallelic, Monoallelic)
        
        Returns:
            True/False for whether the gene and inheritance are consistent with
            the copy number and mechanism.
        """
        
        copy_number_dict = {"DEL": {"0", "1"}, "DUP": {"3"}}
        mech_dict = {"0": {"Uncertain", "Loss of function", \
            "Dominant negative"}, "1": {"Uncertain", "Loss of function", \
            "Dominant negative"}, "3": {"Uncertain", "Increased gene dosage"}}
        
        chrom = "all"
        copy_number = copy_number_dict[self.variant.child.genotype]
        cnv_mech = mech_dict[self.variant.child.info["CNS"]]
        
        if inh == "Biallelic":
            copy_number = {"0"}
        elif inh == "Monoallelic":
            pass
        elif inh == "X-linked dominant":
            chrom = "X"
        elif inh == "Hemizygous":
            chrom = "X"
            if self.variant.child.is_female():
                copy_number = {"3"}
        else:
            # exclude other inheritance modes (Mosaic etc) with impossible
            # criteria
            copy_number = {"XXXX"}
         
        return (chrom =="all" or
            (chrom == "X" and self.variant.get_chrom() == "X")) and \
            self.variant.child.info["CNS"] in copy_number and \
            len(self.known_genes[gene]["inh"][inh] & cnv_mech) > 0
    
    def passes_intragenic_dup(self, gene, inh):
        """ checks if the CNV is an intragenic dup (in an appropriate gene)
        """
        
        allowed_inh = set(["Monoallelic", "X-linked dominant"])
        allowed_mech = set(["Loss of Function"])
        
        # only allow duplications in dominant genes
        if self.variant.child.genotype != "DUP" or inh not in allowed_inh:
            return False
        
        # only allow genes with loss-of-function mechanisms
        if len(self.known_genes[gene]["inh"][inh] & allowed_mech) < 1:
            return False
        
        # find the CNV start and end, as well as the gene start and end
        cnv_start, cnv_end = self.variant.get_range()
        
        gene_start = self.known_genes[gene]["start"]
        gene_end = self.known_genes[gene]["end"]
        
        # check if any part of the gene is outside the CNV boundaries
        return gene_start < cnv_start or gene_end > cnv_end
    
    def check_cnv_region_overlap(self, cnv_regions):
        """ finds CNVs that overlap DECIPHER syndrome regions
        """
        
        chrom = self.variant.child.get_chrom()
        start, end = self.variant.child.get_range()
        copy_number = int(self.variant.child.info["CNS"])
        
        for region_key in cnv_regions:
            region_chrom = region_key[0]
            region_start = int(region_key[1])
            region_end = int(region_key[2])
            region_copy_number = int(cnv_regions[region_key])
            
            if region_chrom != chrom or copy_number != region_copy_number:
                continue
            
            if start <= region_end and end >= region_start and \
                    self.has_enough_overlap(start, end, region_start, region_end):
                return True
        
        return False
    
    def has_enough_overlap(self, start, end, region_start, region_end):
        """ finds if a CNV and another chrom region share sufficient overlap
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
        return forward > 0.01 and reverse > 0.01
        
