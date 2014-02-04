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
    
    def get_candidiate_variants(self):
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
        
        self.check_compound_hets(self.compound_hets)
        
        possibilities = self.candidates + self.compound_hets
        possibilities = self.exclude_duplicates(possibilities)
        
        return possibilities
    
    def exclude_duplicates(self, variants):
        """ rejig variants included under multiple inheritance mechanisms
        
        Args:
            variants: list of candidate variants
        
        Returns:
            list of [variant, check_type, inheritance], with duplicates 
            excluded, and originals modified to show both mechanisms
        """
        
        unique_vars = {}
        for variant in variants:
            key = variant[0].child.get_key()
            if key not in unique_vars:
                unique_vars[key] = variant
            else:
                result = variant[1]
                inh = variant[2]
                
                # append the check type and inheritance type to the first
                # instance of the variant
                if result not in unique_vars[key][1]:
                    unique_vars[key][1] += "," + result
                if inh not in unique_vars[key][2]:
                    unique_vars[key][2] += "," +  inh
        
        return unique_vars.values()
    
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
        
        if check == "compound_het" or check == "hemizygous":
            self.compound_hets.append([variant, check, inheritance])
        elif check == "single_variant":
            self.candidates.append([variant, check, inheritance])
    
    def examine_variant(self, variant, inheritance):
        """ examines a single variant for whether or not to report it
        """
        
        if not self.trio.has_parents():
            return self.check_variant_without_parents(inheritance)
        
        if variant.is_cnv():
            cnv_checker = CNVInheritance(variant, self.trio, self.known_genes)
            check = cnv_checker.check_single_inheritance()
            self.log_string = cnv_checker.log_string
            return check
        
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
            self.compound_hets = []
            return
        
        variant_check = variants[0][2]
        # compound hets on the X chromosome only occur in a female child
        if variant_check == "hemizygous" and self.trio.child.is_male():
            self.compound_hets = []
            return
        
        # check for proband without parental genotypes
        if not self.trio.has_parents():
            return self.compound_hets
        
        use_variants = [False] * len(variants)
        # NB, this is a kind of wasteful way to do this, but since so few
        # vars come through here per sample, it's OK.
        # NB2, de novos plus inh get reported, even though they are ambiguously
        # compound
        for first in range(len(variants)):
            for second in range(len(variants)):
                if first == second:
                    continue
                
                # now we have two different variants in the same gene
                self.set_trio_genotypes(variants[first][0])
                mom_1 = self.mom
                dad_1 = self.dad
                self.set_trio_genotypes(variants[second][0])
                mom_2 = self.mom
                dad_2 = self.dad
                
                # compound hets on the X chromosome occur when the father has a
                # nonref genotype and affected (both might be ref if one is de 
                    # novo, captured later)
                if variant_check == "hemizygous" and \
                   (dad_1.is_hom_alt() or dad_2.is_hom_alt()):
                    if not self.father_affected:
                        continue
                
                # two mutations hitting same gene
                if (mom_1.is_not_ref() and mom_2.is_not_ref() and \
                    dad_1.is_not_ref() and dad_2.is_not_ref()):
                    # if both variants are 1/1/1 hets, we want both parents to 
                    # be affected
                    if self.mother_affected and self.father_affected:
                        use_variants[first] = True
                        use_variants[second] = True
                elif ((mom_1.is_hom_ref() and dad_1.is_hom_ref()) or \
                    (mom_2.is_hom_ref() and dad_2.is_hom_ref())):
                    # one is de novo, so they both definitely get reported
                    use_variants[first] = True
                    use_variants[second] = True
                elif not ((mom_1.is_hom_ref() and mom_2.is_hom_ref()) or \
                    (dad_1.is_hom_ref() and dad_2.is_hom_ref())):
                    use_variants[first] = True
                    use_variants[second] = True
        
        # now add the variants that could be compound hets to a new list
        self.compound_hets = []
        for position in range(len(use_variants)):
            if use_variants[position]:
                # make sure that the check type is compound_het, even for 
                # hemizygous variants
                variants[position][1] = "compound_het"
                self.compound_hets.append(variants[position])
                logging.debug(self.trio.child.get_ID() + " position " + \
                    variants[position][0].get_position() + " is compound het")


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
        
        self.inheritance_modes = set(["X-linked dominant", "Hemizygous", \
            "Monoallelic", "X-linked over-dominance"])
        
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
            report = "hemizygous"
        else:
            raise ValueError("unknown gene inheritance status: " + str(inheritance))
        
        if self.mom.is_hom_ref() and self.dad.is_hom_ref():
            self.log_string = "female x chrom de novo"
            return report
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
        """
        """
        
        self.variant = variant
        self.trio = trio
        self.known_genes = known_genes
        
        self.gene = self.variant.child.gene
    
    def check_single_inheritance(self):
        """ checks if a CNV could be causal by itself
        """
        
        passes = False
        if self.passes_nonddg2p_filter():
            passes = True
        if self.known_genes is not None and self.passes_ddg2p_filter():
            passes = True
            
        if passes:
            return "single_variant"
        else:
            return "nothing"
    
    def check_compound_inheritance(self):
        """ checks if a CNV could contribute to a compound het
        """
        
        passes = False
        if self.passes_nonddg2p_filter():
            passes = True
        if self.known_genes is not None and self.passes_ddg2p_filter():
            passes = True
            
        if passes:
            return "single_variant"
        else:
            return "nothing"
    
    def passes_gene_inheritance(self, gene, inh):
        """ create CNV filter values dependent on DDG2P inheritance mode
        """
        
        # set some default parameters
        chrom = "all"
        # copy_number = (["0", "1", "3"])
        # mechanisms = set(["Uncertain", "Loss of function", \
        #     "Dominant negative", "Increased gene dosage"])
        copy_number = {"DEL": {"0", "1"}, "DUP": {"3"}}
        
        mechanisms = {"0": {"Uncertain", "Loss of function", \
            "Dominant negative"}, "1": {"Uncertain", "Loss of function", \
            "Dominant negative"}, "3": {"Uncertain", "Increased gene dosage"}}
        
        if "Biallelic" == inh:
            copy_number = set(["0"])
            cnv_mech = mechanisms[self.variant.child.info["CNS"]]
        elif "Monoallelic" == inh:
            copy_number = copy_number[self.variant.child.genotype]
            cnv_mech = mechanisms[self.variant.child.info["CNS"]]
        elif "X-linked dominant" == inh:
            chrom = "X"
            copy_number = copy_number[self.variant.child.genotype]
            cnv_mech = mechanisms[self.variant.child.info["CNS"]]
        elif "Hemizygous" == inh and self.variant.child.is_male():
            chrom = "X"
            copy_number = copy_number[self.variant.child.genotype]
            cnv_mech = mechanisms[self.variant.child.info["CNS"]]
        elif "Hemizygous" == inh and self.variant.child.is_female():
            chrom = "X"
            copy_number = set(["3"])
            cnv_mech = mechanisms[self.variant.child.info["CNS"]]
        else:
            # other inheritance modes of "Mosaic", or "Digenic" can be ignored
            # by using impossible criteria
            chrom = "GGG"
            copy_number = set(["999"])
            cnv_mech = set(["Nothing"])
        
        # cnv_encompasses = True
        # if self.variant.child.genotype == "DUP":
        #     if "Increased gene dosage" in self.known_genes[gene]["inheritance"][inh]:
        #         if int(self.variant.get_position()) >= int(self.known_genes[gene]["start"]) or \
        #             int(self.variant.info["END"]) <= int(self.known_genes[gene]["stop"]):
        #             cnv_encompasses = False
        
        print(self.variant, str(chrom), "==", self.variant.get_chrom(), \
            str(copy_number), "==", self.variant.child.info["CNS"], \
            str(cnv_mech), "==", self.known_genes[gene]["inheritance"][inh])
                
        return (chrom =="all" or 
            (chrom == "X" and self.variant.get_chrom() == "X")) and \
            self.variant.child.info["CNS"] in copy_number and \
            len(self.known_genes[gene]["inheritance"][inh] & cnv_mech) > 0 #and \
            # cnv_encompasses
    
    def passes_nonddg2p_filter(self):
        """ checks if a CNV passes the non DDG2P criteria
        """
        
        inh = self.variant.child.format["INHERITANCE"]
        
        if inh == "deNovo" or \
            (inh == "paternal" and self.trio.father.is_affected()) or \
            (inh == "maternal" and self.trio.mother.is_affected()) or \
            (inh == "biparental" and \
            (self.trio.father.is_affected() or self.trio.mother.is_affected())):
            if self.variant.child.genotype == "DEL":
                if float(self.variant.child.info["SVLEN"]) >= 100000:
                    self.log_string = "non-DDG2P, DEL CNV, inh: " + inh
                    return True
            elif self.variant.child.genotype == "DUP":
                if float(self.variant.child.info["SVLEN"]) >= 250000:
                    self.log_string = "non-DDG2P, DUP CNV, inh: " + inh
                    return True
        elif inh not in ["deNovo", "paternal", "maternal", "biparental"]:
            if float(self.variant.child.info["SVLEN"]) >= 500000:
                self.log_string = "non-DDG2P, DUP CNV, unknown inh: " + inh
                return True
        
        # print(self.trio.child.get_ID(), self.variant, \
        #     "nonreported non-DDG2P " + self.variant.child.genotype, \
        #     "CNV, inh=" + inh, "mom_aff=" + str(self.trio.mother.is_affected()),\
        #     str(self.trio.mother.get_affected_status()), \
        #     "len=" + self.variant.child.info["SVLEN"])
        
        self.log_string = "nonreported non-DDG2P " + \
            self.variant.child.genotype + " CNV, inh:" + inh
        return False
    
    def passes_ddg2p_filter(self):
        """ checks if a CNV passes the DDG2P CNV criteria
        """
        
        if self.gene == None:
            return False
        elif "," in self.gene:
            genes = self.gene.split(",")
        else:
            genes  = [self.gene]
        
        self.log_string = "non-reported CNV"
        for gene in genes:
            gene_passes = True
            if gene in self.known_genes:
                self.log_string = "non-reported DDG2P CNV"
                gene_type = self.known_genes[gene]["confirmed_status"]
                
                inh_passes = []
                for inh in self.known_genes[gene]["inheritance"]:
                
                    if "Both DD and IF" in gene_type:
                        self.log_string = "Both DD and IF DDG2P gene"
                        inh_passes.append(True)
                    elif "Confirmed DD Gene" in gene_type or "Probable DD gene" in gene_type:
                        inh_passes.append(self.passes_gene_inheritance(gene, inh))
                
                if not any(inh_passes):
                    gene_passes = False
                
            else:
                gene_passes = False
            
            if gene_passes:
                self.log_string = "DDG2P CNV"
                return True
        
        return False
    
    
    


