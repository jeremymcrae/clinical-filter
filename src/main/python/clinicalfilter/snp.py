""" simple class for handling biallelic variants
"""

class snp(object):
    """ a class to take a genotype coded as 0, 1 or 2, and be able to perform simple functions, like
    reporting whether it is heterozygous, homozygous, or neither, depending on whether the variant 
    is on the X chromosome, and if so, whether the individual is male or female.
    """
    
    def __init__(self, genotype, inheritance_type, gender):
        
        self.genotype = genotype
        self.inheritance_type = inheritance_type
        self.gender = gender
        
        self.ref_allele = 'allele_a'
        self.alt_allele = 'allele_b'
        
        # define some codes used in ped files to identify male and female sexes
        self.male_codes = set(["1", "m", "M", "male"])
        self.female_codes = set(["2", "f", "F", "female"])
        
        self.set_reference_genotypes()
        self.convert_genotype_code_to_alleles()
        
    def convert_genotype_code_to_alleles(self):
        """ converts a 0/1/2 genotype code to a set of alleles, depending on the chromosome
        """
        
        if self.inheritance_type == "autosomal":
            self.convert_autosomal_genotype_code_to_alleles()
        elif self.inheritance_type == "XChrMale" or self.inheritance_type == "XChrFemale":
            self.convert_allosomal_genotype_code_to_alleles()
    
    def set_reference_genotypes(self):
        """ sets reference genotypes for homozygotes and heterozygotes, for checking against
        """
        
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale':
            self.hom_ref = set([self.ref_allele, self.ref_allele])
            self.het = set([self.ref_allele, self.alt_allele])
            self.hom_alt = set([self.alt_allele, self.alt_allele])
        elif self.inheritance_type == 'XChrMale':
            self.hom_alt = set([self.alt_allele])
            self.hom_ref = set([self.ref_allele])
    
    def is_het(self):
        """ returns whether a variant is heterozygous
        """
        
        is_het = False
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale':
            if self.alleles == self.het:
                is_het = True
        elif self.inheritance_type == 'XChrMale':
            pass
        
        return is_het
    
    def is_hom_alt(self):
        """ returns whether a genotype is homozygous for the alternate allele
        """
        
        is_hom_alt = False
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale' or self.inheritance_type == 'XChrMale':
            if self.alleles == self.hom_alt:
                is_hom_alt = True
        
        return is_hom_alt
    
    def is_hom_ref(self):
        """ returns whether a variant is homozygous for the reference allele
        """
        
        is_hom_ref = False
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale' or self.inheritance_type == 'XChrMale':
            if self.alleles == self.hom_ref:
                is_hom_ref = True
        
        return is_hom_ref
    
    def is_not_ref(self):
        """ returns whether a variant is not the homozygous for the reference allele
        """
        
        is_not_hom_ref = False
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale' or self.inheritance_type == 'XChrMale':
            if self.alleles != self.hom_ref:
                is_not_hom_ref = True
        
        return is_not_hom_ref
    
    def is_not_alt(self):
        """ returns whether a variant is not homozygous for the alternate allele
        """
        
        is_not_hom_alt = False
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale' or self.inheritance_type == 'XChrMale':
            if self.alleles != self.hom_alt:
                is_not_hom_alt = True
        
        return is_not_hom_alt
    
    def convert_autosomal_genotype_code_to_alleles(self):
        """converts a genotype code to a set of alleles
        
        returns:
            an unsorted set of allele codes
        """
        
        genotype = str(self.genotype)
        
        if genotype == '0':
            self.alleles = set([self.ref_allele, self.ref_allele])
        elif genotype == '1':
            self.alleles = set([self.ref_allele, self.alt_allele])
        elif genotype == "2":
            self.alleles = set([self.alt_allele, self.alt_allele])
        else:
            raise ValueError("genotype code '" + genotype + "' is not recognised")
        
    def convert_allosomal_genotype_code_to_alleles(self):
        """converts a genotype on the x-chromosome into the set of alleles
        
        returns:
            an unsorted set of allele codes
        """
        
        genotype = str(self.genotype)
        
        if self.gender in self.male_codes:
            if genotype == "0":
                self.alleles = set([self.ref_allele])
            elif genotype == "2":
                self.alleles = set([self.alt_allele])
            elif genotype == "1":
                self.alleles = set([self.ref_allele, self.alt_allele])
                raise ValueError("males shouldn't be heterozygous on a sex chromosome")
            else:
                raise ValueError("genotype code '" + str(genotype) + "' is not recognised")
        elif self.gender in self.female_codes:
            self.convert_autosomal_genotype_code_to_alleles()
        else:
            raise ValueError("Unknown gender: " + self.gender)


