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

from clinicalfilter.variant.variant import Variant

class SNV(Variant):
    """ a class to take a SNV genotype for an individual, and be able to perform
    simple functions, like reporting whether it is heterozygous, homozygous, or
    neither, depending on whether the variant is on the X chromosome, and if so,
    whether the individual is male or female.
    """
    
    debug_chrom = None
    debug_pos = None
    
    @classmethod
    def set_debug(cls_obj, chrom, pos):
        cls_obj.debug_chrom = chrom
        cls_obj.debug_pos = pos
        
        if cls_obj.debug_chrom is not None:
            cls_obj.passes_filters = cls_obj.passes_filters_with_debug
    
    def is_cnv(self):
        """ checks whether the variant is for a CNV
        """
        
        return False
    
    def get_key(self):
        """ return a tuple to identify the variant
        """
        
        return (self.get_chrom(), self.get_position())
    
    def set_genotype(self):
        """ sets the genotype of the variant using the format entry
        """
        
        if self.format is not None:
            self.genotype = self.convert_genotype(self.format["GT"])
        else:
            raise ValueError("cannot find a genotype")
        
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
        elif self.inheritance_type in ["XChrMale", "XChrFemale"]:
            self.convert_allosomal_genotype_code_to_alleles()
    
    def set_reference_genotypes(self):
        """ sets reference genotypes for homozygotes and heterozygotes
        """
        
        if self.inheritance_type in ["autosomal", "XChrFemale"]:
            self.hom_ref = set([self.ref_allele, self.ref_allele])
            self.het = set([self.ref_allele, self.alt_alleles])
            self.hom_alt = set([self.alt_alleles, self.alt_alleles])
        elif self.inheritance_type == "XChrMale":
            self.hom_alt = set([self.alt_alleles])
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
            self.alleles = set([self.ref_allele, self.alt_alleles])
        elif genotype == "2":
            self.alleles = set([self.alt_alleles, self.alt_alleles])
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
                self.alleles = set([self.alt_alleles])
            elif genotype == "1":
                raise ValueError("heterozygous X-chromomosome male")
            else:
                raise ValueError("unknown genotype '" + str(genotype))
        else:
            self.convert_autosomal_genotype_code_to_alleles()
    
    def passes_filters(self):
        """Checks whether a VCF record passes user defined criteria.
            
        Returns:
            boolean value for whether the variant passes the filters
        """
        
        pass_value, key = self.check_filters()
        
        return pass_value
    
    def passes_filters_with_debug(self):
        """Checks whether a VCF record passes user defined criteria.
        
        This method replaces passes_filters() when we specify a chromosome and
        position for debugging the filtering.
            
        Returns:
            boolean value for whether the variant passes the filters
        """
        
        pass_value, key = self.check_filters()
        
        if pass_value == False and self.get_position() == self.debug_pos:
            
            if key == "MAF":
                value = self.find_max_allele_frequency()
            elif key == "consequence":
                value = self.consequence
            elif key == "FILTER":
                value = self.filter
            elif key == "HGNC":
                value = "not in a known gene"
            
            print("failed {0}: {1}".format(key, value))
        
        return pass_value
    
    def check_filters(self):
        """Checks whether a VCF record passes user defined criteria.
        
        Returns:
            tuple of (True/False for whether the variant passes the filters, and
                string for the last checked filter)
        """
        
        # exclude variants without functional consequences
        if not self.is_lof() and not self.is_missense():
            return (False, "consequence")
        
        # exclude variants with high minor allele frequencies in any population
        max_maf = self.find_max_allele_frequency()
        if max_maf is not None and max_maf > 0.005:
            return (False, "MAF")
        
        # exclude variants outside genes known to be involved in genetic
        # disorders, unless there isn't any such set of genes available
        genes = [ l for sublist in self.get_genes() for l in sublist ]
        if self.known_genes is not None and len(set(genes) & set(self.known_genes)) == 0:
            return (False, "HGNC")
        
        # exclude variants without PASS values, except where the fail reason is
        # low_VQSLOD and the variant has been detected by denovogear
        if self.filter not in ["PASS", ".", "LOW_VQSLOD"]:
            if ("DENOVO-SNP" not in self.info and "DENOVO-INDEL" not in self.info):
                return (False, "FILTER")
        
        return (True, "passed all")
