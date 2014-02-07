""" A class for checking whether the genotypes of a trio for a variant or 
variants in a gene fit an inheritance model specific to the gene and 
chromosome that the variant/s are in.
"""

import logging


class Inheritance(object):
    """figure out whether trio genotypes fit mendelian inheritance models
    """
    
    def __init__(self, variants, trio, known_genes, gene_inheritance=None):
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
            variants: list of variants in a gene
            trio: a Family object for the family
            gene_inheritance: inheritance based on a database of genes known to
                be involved in disorders. If no gene list is provided, default
                to including every possible inheritance type.
        """
        
        self.variants = variants
        self.trio = trio
        self.known_genes = known_genes
        
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
            
            # check against every inheritance mode for the gene
            for inheritance in self.inheritance_modes & self.gene_inheritance:
                check = self.examine_variant(variant, inheritance)
                self.add_variant_to_appropriate_list(variant, check, inheritance)
            
            logging.debug(self.trio.child.get_ID() + " position " + \
                variant.get_position() + " " + self.log_string)
        
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
        """
        
        if check == "compound_het":
            self.compound_hets.append((variant, check, inheritance))
        elif check == "single_variant":
            self.candidates.append((variant, check, inheritance))
    
    def examine_variant(self, variant, inheritance):
        """ examines a single variant for whether or not to report it
        """
        
        if variant.is_cnv():
            cnv_checker = CNVInheritance(variant, self.trio, self.known_genes)
            check = cnv_checker.check_single_inheritance()
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
    
    def __init__(self, variants, trio, known_genes, gene_inheritance=None):
        
        super(Autosomal, self).__init__(variants, trio, known_genes, gene_inheritance)
        
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
            self.log_string = "de novo"
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
        
        self.log_string = "non-causal homozygous variant"
        return "nothing"


class Allosomal(Inheritance):
    
    def __init__(self, variants, trio, known_genes, gene_inheritance=None):
        
        super(Allosomal, self).__init__(variants, trio, known_genes, gene_inheritance)
        
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
        
        elif self.trio.child.is_female():
            if self.dad.is_hom_ref() or self.mom.is_hom_ref():
                self.log_string = "female child hom alt and father hom ref, which is non-mendelian"
                return "nothing"
            elif (self.mom.is_het() or \
                 (self.mom.is_hom_alt() and self.mother_affected)) and \
                 (self.dad.is_hom_alt() and self.father_affected):
                self.log_string = "testing"
                return "single_variant"
        
        self.log_string = "variant not compatible with being causal"
        return "nothing"
   

class CNVInheritance(object):
    
    def __init__(self, variant, trio, known_genes):
        """ intialise the class
        
        Args: 
            variant: CNV to check for inheritance
            trio: family trio object
            known_genes: dictionary of known genes, currently the DDG2P set
        """
        
        self.variant = variant
        self.trio = trio
        self.known_genes = known_genes
    
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
        # elif self.check_compound_inheritance():
        #     return "compound_het"
        
        if self.check_compound_inheritance():
            self.log_string = "possible compound het CNV"
            return "compound_het"
        
        return "nothing"
    
    def check_compound_inheritance(self):
        """ checks if a CNV could contribute to a compound het
        """
        
        # we don't want CNVs that don't have copy number of 1 or 3, since 
        # copy number = 1 or 3 are the only ones that could operate as compound
        # hets (other copy numbers such as 0 are implicitly dominant)
        if self.variant.child.info["CNS"] not in {"1", "3"}:
            return False
            # return True
        
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
                if "Biallelic" in self.known_genes[gene] or \
                    (self.variant.get_chrom() == "X" and \
                    "Hemizygous" in self.known_genes[gene] and \
                    self.trio.child.is_female() and \
                    self.variant.child.info["CNS"] == "1"):
                   return True
        
        # return True
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
        if inh in ["deNovo", "paternal", "maternal", "biparental", "inheritedDuo"]:
            if geno == "DEL":
                min_len = 100000
            elif geno == "DUP":
                min_len = 250000
        else:
            min_len = 500000
        
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
            gene_type = self.known_genes[gene]["confirmed_status"]
            
            inh_passes = []
            if "Both DD and IF" in gene_type:
                self.log_string = "Both DD and IF DDG2P gene"
                inh_passes.append(True)
            elif {"Confirmed DD Gene", "Probable DD gene"} & gene_type == 0:
                continue
            
            for inh in self.known_genes[gene]["inheritance"]:
                inh_passes.append(self.passes_gene_inheritance(gene, inh))
            
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
        
        # cnv_encompasses = True
        # if self.variant.child.genotype == "DUP":
        #     if "Increased gene dosage" in self.known_genes[gene]["inheritance"][inh]:
        #         if int(self.variant.get_position()) >= int(self.known_genes[gene]["start"]) or \
        #             int(self.variant.child.info["END"]) <= int(self.known_genes[gene]["stop"]):
        #             cnv_encompasses = False
         
        return (chrom =="all" or 
            (chrom == "X" and self.variant.get_chrom() == "X")) and \
            self.variant.child.info["CNS"] in copy_number and \
            len(self.known_genes[gene]["inheritance"][inh] & cnv_mech) > 0 # and \
            # cnv_encompasses
    
    


