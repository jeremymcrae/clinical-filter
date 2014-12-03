""" class for holding single nucleotide variant data for a single individual
"""

from clinicalfilter.vcf_info import VcfInfo
from clinicalfilter.variant import Variant

class SNV(Variant, VcfInfo):
    """ a class to take a SNV genotype for an individual, and be able to perform 
    simple functions, like reporting whether it is heterozygous, homozygous, or 
    neither, depending on whether the variant is on the X chromosome, and if so,
    whether the individual is male or female.
    """
    
    def get_key(self):
        """ return a tuple to identify the variant
        """
        
        return (self.get_chrom(), self.get_position())
    
    def set_genotype(self):
        """ sets the genotype of the variant using the format entry
        """
        
        if hasattr(self, "format"):
            self.genotype = self.convert_genotype(self.format["GT"])
        else:
            raise ValueError("cannot find a genotype")
        
        self.set_reference_genotypes()
        self.convert_genotype_code_to_alleles()
    
    def set_default_genotype(self):
        """ for variants lacking genotypes, set a default genotype
        """
        
        self.genotype = 0
        
        self.set_reference_genotypes()
        self.convert_genotype_code_to_alleles()
    
    def convert_genotype(self, genotype):
        """Maps genotypes from two character format to single character.
        
        Args:
            genotype: genotype in two character format. eg "0/0"
        
        Returns:
            Count of non-reference alleles
        """
        
        if len(genotype) == 1:
            raise ValueError("genotype is only a single character")
        
        # split the genotype field (allow for phased genotypes)
        try:
            allele_1, allele_2 = genotype.split("/")
        except ValueError:
            allele_1, allele_2 = genotype.split("|")
        
        assert self.is_number(allele_1)
        assert self.is_number(allele_2)
        
        # if the two alleles are different, return 1, which roughly equates
        # to heterozygous. Strictly this isn't quite true, since some variants
        # might have both alleles as non-reference, but different from each
        # other. The cases where this occurs all occur for indels, and appear to
        # be poorly called variants, where it is likely that one of the alleles
        # is actually for the reference. 
        if allele_1 != allele_2:
            return 1
        elif allele_1 == "0" and allele_2 == "0":
            return 0
        
        return 2
    
    def convert_genotype_code_to_alleles(self):
        """ converts a genotype to a set of alleles
        """
        
        if self.inheritance_type == "autosomal":
            self.convert_autosomal_genotype_code_to_alleles()
        elif self.inheritance_type == "XChrMale" or self.inheritance_type == "XChrFemale":
            self.convert_allosomal_genotype_code_to_alleles()
    
    def set_reference_genotypes(self):
        """ sets reference genotypes for homozygotess and heterozygotes
        """
        
        if self.inheritance_type == "autosomal" or self.inheritance_type == "XChrFemale":
            self.hom_ref = set([self.ref_allele, self.ref_allele])
            self.het = set([self.ref_allele, self.alt_allele])
            self.hom_alt = set([self.alt_allele, self.alt_allele])
        elif self.inheritance_type == "XChrMale":
            self.hom_alt = set([self.alt_allele])
            self.het = set([])
            self.hom_ref = set([self.ref_allele])
        else:
            raise ValueError("unknown inheritance type:", self.inheritance_type)
    
    def is_het(self):
        """ returns whether a variant is heterozygous
        """
        
        return self.alleles == self.het
    
    def is_hom_alt(self):
        """ returns whether a genotype is homozygous for the alternate allele
        """
        
        return self.alleles == self.hom_alt
    
    def is_hom_ref(self):
        """ returns whether a variant is homozygous for the reference allele
        """
        
        return self.alleles == self.hom_ref
    
    def is_not_ref(self):
        """ returns whether a variant is not homozygous for the reference allele
        """
        
        return self.alleles != self.hom_ref
    
    def is_not_alt(self):
        """ returns whether a variant is not homozygous for the alternate allele
        """
        
        return self.alleles != self.hom_alt
    
    def convert_autosomal_genotype_code_to_alleles(self):
        """converts a genotype code to a set of alleles
        """
        
        genotype = str(self.genotype)
        
        if genotype == "0":
            self.alleles = set([self.ref_allele, self.ref_allele])
        elif genotype == "1":
            self.alleles = set([self.ref_allele, self.alt_allele])
        elif genotype == "2":
            self.alleles = set([self.alt_allele, self.alt_allele])
        else:
            raise ValueError("unknown genotype '" + str(genotype))
        
    def convert_allosomal_genotype_code_to_alleles(self):
        """converts a genotype on the x-chromosome into the set of alleles
        """
        
        genotype = str(self.genotype)
        
        if self.is_male():
            if genotype == "0":
                self.alleles = set([self.ref_allele])
            elif genotype == "2":
                self.alleles = set([self.alt_allele])
            elif genotype == "1":
                raise ValueError("heterozygous X-chromomosome male")
            else:
                raise ValueError("unknown genotype '" + str(genotype))
        elif self.is_female():
            self.convert_autosomal_genotype_code_to_alleles()
        else:
            raise ValueError("Unknown gender: " + self.gender)
