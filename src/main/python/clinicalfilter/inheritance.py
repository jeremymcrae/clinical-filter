""" A class for checking whether the genotypes of a trio for a variant or 
variants in a gene fit an inheritance model specific to the gene and 
chromosome that the variant/s are in.
"""

import logging


class Inheritance(object):
    """figure out whether trio genotypes fit mendelian inheritance models
    """
    
    def __init__(self, variants, trio, gene_inheritance=None):
        """ intialise the class with the variants and trio information
        
        We have an affected child, with two parents who may or may not be affected.
        we want a model where either some alleles are transmitted to the child, which might 
        contribute to mendelian inheritance of a disorder, or the child has a de novo mutation 
        that might contribute to the child's disorder.
        Genotypes are encoded as 0, 1, or 2, which indicate the number of non reference alleles 
        (where reference appears to be arbitrarily defined as the major allele in a population, 
        using the 1000 Genomes populations in this instance).
        The genotypes on the sex chromosomes are still defined as 0, 1, and 2. Males with a
        non-reference genotype on the x-chromosome are coded as 2, and females with a
        non-reference genotype on the Y-chromosome are coded as 2. The genotypes that are 
        processed through here should be for rare variants, as defined as having a minor allele 
        frequency less than 0.01 in populations (typically 1000 Genomes populations, as well as 
        any internal reference populations).
        
        Also, account for X inactivation (X chromosome mosaicism in females) and non-mendelian 
        inheritance. Non-mendelian inheritance can either occur from de novo mutation, or from 
        genotype miscalls. De novo mutations are more likely to be reported in genome positions 
        that are not already variable (due to the fact that only 1% of the genome has been 
        reported to contain SNPs.) Apparent de novo mutations at previously identified variable 
        sites are therefore more likely to be miscalls.
        
        also, for some genes we know whether the gene is monoallelic, or biallelic, or other 
        possibilities. This also inform whether genotypes for a trio could be causal.
        
        Args:
            variants: list of variants in a gene
            trio: a pedTrio object for the family
            gene_inheritance: inheritance based on a database of genes known to be involved in 
                disorders. If no gene list is provided, default to including every possible 
                inheritance type.
        """
        
        self.variants = variants
        self.trio = trio
        
        if self.trio.father is not None:
            self.father_affected = self.trio.father.is_affected()
        else:
            self.father_affected = None
            
        if self.trio.mother is not None:
            self.mother_affected = self.trio.mother.is_affected()
        else:
            self.mother_affected = None
        
        self.chrom_inheritance = self.variants[0].get_inheritance_type()
        self.gene_inheritance = gene_inheritance
        
        # here are the complete set of inheritance modes defined in the knwon gene database
        if self.gene_inheritance is None:
            self.gene_inheritance = set(["Biallelic", "Both", "Digenic", "Hemizygous", "Imprinted", \
                                    "Mitochondrial", "Monoallelic", "Mosaic", "Uncertain", \
                                    "X-linked dominant", "X-linked over-dominance"])
        
        # if a gene has an inheritance mode of "Both", make sure we will process both mono
        # and biallelic inheritance, but remove "Both"
        if "Both" in self.gene_inheritance:
            self.gene_inheritance.add("Biallelic")
            self.gene_inheritance.add("Monoallelic")
            self.gene_inheritance.remove("Both")
    
    def get_candidiate_variants(self):
        """ screen for variants that might contribute to an affected childs disorder
        """
        
        if self.check_inheritance_mode_matches_gene_mode() == False:
            return []
         
        if self.trio.father is None and self.trio.mother is None:
            possibilities = self.find_variants_without_parents()
        else:
            possibilities = self.find_variants_for_trio()
        
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
        
        variant_keys = set()
        variants_to_keep = []
        for variant in variants:
            key = variant[0].child.get_key()
            if key not in variant_keys:
                variants_to_keep.append(variant)
                variant_keys.add(key)
            else:
                for pos in range(len(variants_to_keep)):
                    temp = variants_to_keep[pos]
                    if temp[0].child.get_key() == key:
                        temp[1] += "," + variant[1]
                        temp[2] += "," +  variant[2]
                        variants_to_keep[pos] = temp
        
        return variants_to_keep
    
    def check_inheritance_mode_matches_gene_mode(self):
        """ make sure that the mode of inheritance for the gene makes sense 
        for the mode of inheritance for the chromosome - there's no point 
        looking through hemizygous genes on an autosomal chromosome, or 
        biallelic genes on an X chromosome.
        """
        
        # at this point make sure we only deal with genes that have at least 
        # one correct modes of inheritance for the given chromosome type
        if len(self.inheritance_modes & self.gene_inheritance) > 0:
            return True
        
        return False
    
    def set_trio_genotypes(self, variant):
        """ sets the genotypes for the trio as class objects
        """
        
        # allow for children without parents
        if self.trio.father is None and self.trio.mother is None:
            self.child = variant.child
            self.dad = None
            self.mom = None
        else:
            self.child = variant.child
            self.mom = variant.mother
            self.dad = variant.father
    
    def passes_sift_and_polyphen(self, variant):
        """ returns whether a variant has sift and polyphen annotations for minimal impact
        """
        
        if "SIFT" not in variant.info or "PolyPhen" not in variant.info:
            return True
        elif "tolerated" in variant.info["SIFT"] and "benign" in variant.info["PolyPhen"]:
            return False
        else:
            return True
    
    def add_variant_to_appropriate_list(self, variant, variant_check, inheritance_type):
        """ add processed variants to the appropriate list
        """
        
        if variant_check == "compound_het" or variant_check == "hemizygous":
            self.compound_hets.append([variant, variant_check, inheritance_type])
        elif variant_check == "single_variant":
            self.candidates.append([variant, variant_check, inheritance_type])
    
    def find_variants_for_trio(self):
        """ find variants when we have trio (child, mother, father) information
        """
        
        self.compound_hets = []
        self.candidates = []
        
        for variant in self.variants:
            self.nucleotide_position = variant.get_position()
            
            # set the genotypes for the trio, unless any of the genotypes are unsettable (due to
            # being male X hets etc) in which case, just pass on to the next variant
            try:
                self.set_trio_genotypes(variant)
            except ValueError:
                logging.warning(self.trio.child.get_ID() + " position " + str(self.nucleotide_position) + ": unable to use genotypes")
                continue
            
            # allow for genes having multiple inheritance modes by analysing them separately
            for inheritance in self.inheritance_modes & self.gene_inheritance:
                variant_check = self.examine_variant(variant, inheritance)
                self.add_variant_to_appropriate_list(variant, variant_check, inheritance)
            
            logging.debug(self.trio.child.get_ID() + " position " + str(self.nucleotide_position) + " " + self.log_string)
        
        if len(self.compound_hets) > 0:
            self.compound_hets = self.check_compound_hets(self.compound_hets)
        
        possibilities = self.candidates + self.compound_hets
        return possibilities
    
    def examine_variant(self, variant, inheritance):
        """ examines a single variant for whether or not to report it
        """
        
        if self.child.is_het():
            return self.check_heterozygous_variant(inheritance)
        elif self.child.is_hom_alt():
            if self.check_homozygous_variant(inheritance):
                return "single_variant"
        else:
            self.log_string = "neither homozygous alt nor heterozygous: " + str(self.child) + " with alleles: " + str(self.child) +  "and inheritance type:" + self.chrom_inheritance
            return "nothing"
            
        logging.error(self.trio.child.get_ID() + " position " + str(self.nucleotide_position) + " something has gone wrong here")
    
    def check_compound_hets(self, variants):
        """ checks for compound hets within a gene
        """
        
        variant_check = variants[0][2]
        # compound hets on the X chromosome only occur in a female child
        if variant_check == "hemizygous" and self.trio.child.is_male():
            return []
        
        use_variants = [False] * len(variants)
        # NB, this is a kind of wasteful way to do this, but since so few
        # vars come through here per sample, it's OK.
        # NB2, de novos plus inh get reported, even though they are ambiguously compound
        for first in range(len(variants)):
            for second in range(len(variants)):
                if first == second:
                    continue
                
                # ok, so we are in the same gene, and looking at two different variants
                self.set_trio_genotypes(variants[first][0])
                mom_1 = self.mom
                dad_1 = self.dad
                self.set_trio_genotypes(variants[second][0])
                mom_2 = self.mom
                dad_2 = self.dad
                
                # compound hets on the X chromosome occur when the father has a nonref genotype 
                # and affected (both might be ref if one is de novo, captured later)
                if variant_check == "hemizygous" and (dad_1.is_hom_alt() or dad_2.is_hom_alt()):
                    if not self.father_affected:
                        continue
                
                # two mutations hitting same gene
                if (mom_1.is_not_ref() and mom_2.is_not_ref() and dad_1.is_not_ref() and dad_2.is_not_ref()):
                    # if both variants are 1/1/1 hets, we want both parents to be affected
                    if self.mother_affected and self.father_affected:
                        use_variants[first] = True
                        use_variants[second] = True
                elif ((mom_1.is_hom_ref() and dad_1.is_hom_ref()) or (mom_2.is_hom_ref() and dad_2.is_hom_ref())):
                    # one is de novo, so they both definitely get reported
                    use_variants[first] = True
                    use_variants[second] = True
                elif not ((mom_1.is_hom_ref() and mom_2.is_hom_ref()) or (dad_1.is_hom_ref() and dad_2.is_hom_ref())):
                    use_variants[first] = True
                    use_variants[second] = True
        
        # now add the variants that could be compound hets to a new list
        compound_hets = []
        for position in range(len(use_variants)):
            if use_variants[position]:
                compound_hets.append(variants[position])
                logging.debug(self.trio.child.get_ID() + " position " + str(variants[position][1]) + " is compound het")
        
        return compound_hets


class Autosomal(Inheritance):
    
    def __init__(self, variants, trio, gene_inheritance=None):
        
        super(Autosomal, self).__init__(variants, trio, gene_inheritance)
        
        self.inheritance_modes = set(["Monoallelic", "Biallelic", "Both"])
    
    def find_variants_without_parents(self):
        """ test variants in children where we lack parental genotype information.
        """
        
        self.compound_hets = []
        self.candidates = []
        
        for variant in self.variants:
            self.nucleotide_position = variant.get_position()
            # set the genotypes for the trio, unless any of the genotypes are unsettable (due to
            # being male X hets etc) in which case, just pass on to the next variant
            try:
                self.set_trio_genotypes(variant)
            except ValueError:
                logging.warning(self.trio.child.get_ID() + " position " + str(self.nucleotide_position) + ": unable to use genotypes")
                continue
            
            for inheritance in self.inheritance_modes & self.gene_inheritance:
                if self.child.is_het() and inheritance == "Biallelic":
                    self.add_variant_to_appropriate_list(variant, "compound_het", inheritance)
                elif self.child.is_hom_alt() and inheritance == "Biallelic":
                    self.add_variant_to_appropriate_list(variant, "single_variant", inheritance)
                elif self.child.is_het() and inheritance == "Monoallelic":
                    self.add_variant_to_appropriate_list(variant, "single_variant", inheritance)
        
        # if we have less than two hets in the compound het candidates, then that's not enough to
        # make a compound het
        if len(self.compound_hets) < 2:
            self.compound_hets = [] 
        
        possibilities = self.candidates + self.compound_hets
        
        return possibilities
    
    def check_heterozygous_variant(self, inheritance):
        """ checks, for a het child, if a trio's genotypes could contribute to disease
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
        elif (self.dad.is_not_ref() and self.father_affected) and (self.mom.is_hom_ref() or self.mother_affected) or \
             (self.mom.is_not_ref() and self.mother_affected) and (self.dad.is_hom_ref() or self.father_affected):
            self.log_string = "transmitted from aff, other parent non-carrier or aff"
            return report
        elif "Biallelic" == inheritance and ((self.dad.is_not_alt() or self.father_affected) and (self.mom.is_not_alt() or self.mother_affected)):
            self.log_string = "het-check for recessive genes and unaff parents not homoz"
            return report
        else:
            self.log_string = "typically for trios with non-de novo unaffected parents"
        
    def check_homozygous_variant(self, inheritance):
        """ checks, for a homozygous child, if a trio's genotypes could contribute to disease
        """
        
        if self.dad.is_hom_ref() or self.mom.is_hom_ref():
            #NB: will miss possibility of one inherited copy and one de-novo at same site.
            self.log_string = "child hom alt and parents hom ref, which is non-mendelian"
            return False
        elif "Biallelic" == inheritance:
            if (self.mom.is_het() and self.dad.is_het()):
                self.log_string = "both parents het"
                return True
            elif (((self.mom.is_hom_alt() and self.mother_affected) and (self.dad.is_not_alt() or self.father_affected)) or \
                ((self.dad.is_hom_alt() and self.father_affected) and (self.mom.is_not_alt() or self.mother_affected))):
                self.log_string = "homoz parent aff"
                return True
        elif "Monoallelic" == inheritance:
            # dominant
            if (((self.dad.is_not_ref() and self.father_affected) and (self.mom.is_hom_ref() or self.mother_affected)) or \
               ((self.mom.is_not_ref() and self.mother_affected) and (self.dad.is_hom_ref() or self.father_affected))):
                self.log_string = "transmitted from aff, other parent non-carrier or aff"
                return True
        
        self.log_string = "variant not compatible with being causal"
        return False


class Allosomal(Inheritance):
    
    def __init__(self, variants, trio, gene_inheritance=None):
        
        super(Allosomal, self).__init__(variants, trio, gene_inheritance)
        
        self.inheritance_modes = set(["X-linked dominant", "Hemizygous", \
            "Monoallelic", "X-linked over-dominance"])
        
        # on the X chrom, treat monoallelic and X-linked dominant modes of 
        # inheritance the same
        if "Monoallelic" in self.gene_inheritance:
            self.gene_inheritance.add("X-linked dominant")
            self.gene_inheritance.remove("Monoallelic")
    
    def find_variants_without_parents(self):
        """ test variants in children where we lack parental genotype information.
        """
        
        self.compound_hets = []
        self.candidates = []
        
        for variant in self.variants:
            self.nucleotide_position = variant.get_position()
            # set the genotypes for the trio, unless any of the genotypes are unsettable (due to
            # being male X hets etc) in which case, just pass on to the next variant
            try:
                self.set_trio_genotypes(variant)
            except ValueError:
                logging.warning(self.trio.child.get_ID() + " position " + str(self.nucleotide_position) + ": unable to use genotypes")
                continue
            
            # allow for genes having multiple inheritance modes by analysing them separately
            for inheritance in self.inheritance_modes & self.gene_inheritance:
                if (self.child.is_hom_alt() and inheritance == "X-linked dominant") or \
                   (self.trio.child.is_female() and self.child.is_het() and inheritance == "X-linked dominant"):
                    self.add_variant_to_appropriate_list(variant, "single_variant", inheritance)
                elif (self.child.is_hom_alt() and inheritance == "Hemizygous"):
                    self.add_variant_to_appropriate_list(variant, "single_variant", inheritance)
                elif (self.trio.child.is_female() and self.child.is_het() and inheritance == "Hemizygous"):
                    self.add_variant_to_appropriate_list(variant, "hemizygous", inheritance)
        
        # if we have less than two hets in the compound het candidates, then that's not enough to
        # make a compound het
        if len(self.compound_hets) < 2:
            self.compound_hets = [] 
        
        possibilities = self.candidates + self.compound_hets
        
        return possibilities
    
    def check_heterozygous_variant(self, inheritance):
        """ checks, for a X chrom heterozygote, whether a trios genotypes could contribute
        """
        
        if self.trio.child.is_male():
            self.log_string = "male proband should not be heterozygous on the X chromosome"
            return False
        
        if "X-linked dominant" == inheritance or "Monoallelic" == inheritance:
            report = "single_variant"
        elif "Hemizygous" == inheritance:
            # recessive: should be marked for compound-het screen
            report = "hemizygous"
        elif "X-linked over-dominance" == inheritance:
            report = "X-linked over dominance not currently supported"
            return "nothing"
        else:
            raise ValueError("unknown gene inheritance status: " + str(inheritance))
        
        if self.mom.is_hom_ref() and self.dad.is_hom_ref():
            self.log_string = "female x chrom de novo"
            return report
        elif (self.dad.is_hom_alt() and self.father_affected) and (self.mom.is_hom_ref() or self.mother_affected) or \
           (self.mom.is_not_ref() and self.mother_affected) and (self.dad.is_hom_ref() or self.father_affected):
            self.log_string = "x chrom transmitted from aff, other parent non-carrier or aff"
            return report
        elif inheritance == "X-linked over-dominance" and not self.father_affected and \
             ((self.mom.is_hom_ref() and self.dad.is_hom_alt()) or \
             self.mom.is_het() and self.mother_affected):
            self.log_string = "X-linked inheritance with unaffected hom alt males and females, but affected het females (eg PCDH19)"
            return "single_variant"
        else:
            self.log_string = "variant not compatible with being causal"
            return "nothing"
        
    def check_homozygous_variant(self, inheritance):
        """ checks, for a X chrom homozygote, whether a trios genotypes could contribute
        """
        
        if "X-linked dominant" == inheritance:
            report = "single_variant"
        elif "Hemizygous" == inheritance:
            # recessive: should be marked for compound-het screen
            report = "hemizygous"
        elif inheritance == "X-linked over-dominance":
            self.log_string = "X-linked over dominance not currently supported"
            return False
        else:
            raise ValueError("unknown gene inheritance status: " + str(inheritance))
        
        # treat male sex inheritance differently from female sex inheritance
        if self.trio.child.is_male():
            if self.mom.is_hom_ref():
                self.log_string = "male X chrom de novo"
                return True
            elif (self.mom.is_het() and not self.mother_affected) or (self.mom.is_hom_alt() and self.mother_affected):
                self.log_string = "male X chrom inherited from het mother or hom affected mother"
                return True
        
        elif self.trio.child.is_female():
            if self.dad.is_hom_ref() or self.mom.is_hom_ref():
                self.log_string = "female child hom alt and father hom ref, which is non-mendelian"
                return False
            elif (self.mom.is_het() or (self.mom.is_hom_alt() and self.mother_affected)) and (self.dad.is_hom_alt() and self.father_affected):
                self.log_string = "testing"
                return True
        
        self.log_string = "variant not compatible with being causal"
        return False




