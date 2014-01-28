""" class for holding copy number information for a single individual
"""

from clinicalfilter.vcf_info import VcfInfo
from clinicalfilter.variant import Variant

class CNV(Variant, VcfInfo):
    """  class to take CNV data for an individual, and
    """
    
    def set_genotype(self):
        """ sets the genotype of the variant
        """
        
        self.ref_genotypes = set(["REF"])
        self.alt_genotypes = set(["DEL", "DUP"])
        
        if self.alt_allele == "<DUP>":
            self.genotype = "DUP"
        elif self.alt_allele == "<DEL>":
            self.genotype = "DEL"
        else:
            raise ValueError("unknown CNV allele code")
        
        self.set_range()
    
    def set_default_genotype(self):
        """ set a default genotype for individuals without one
        """
        
        self.genotype = "REF"
        self.set_range()
    
    def set_range(self):
        """ sets the range for the CNV
        """
        
        self.start_position = self.position
        
        if hasattr(self, "info"):
            self.end_position = self.info["END"]
        else:
            self.end_position = str(int(self.start_position) + 10000)
        
        self.range = (self.start_position, self.end_position)
    
    def get_key(self):
        """ return a tuple to identify the variant
        """
        
        try:
            return (self.chrom, self.start_position, self.end_position)
        except AttributeError:
            self.set_range()
        
        return (self.chrom, self.start_position, self.end_position)
    
    def add_gene_from_info(self):
        """ adds a gene to the var using the info. CNVs and SNVs act differently
        """
        
        # sometimes the variant lacks an HGNC field, but does have a HGNC_ALL
        # entry.
        
        if "HGNC_ALL" in self.info:
            self.gene = self.info["HGNC_ALL"]
        elif "HGNC" not in self.info and "NUMBERGENES" in self.info:
            self.gene = None
            if int(self.info["NUMBERGENES"]) > 0:
                self.gene = "."
        # make sure we have gene and consequence keys in the info dict, for 
        # the filter to work with
        elif "HGNC" not in self.info:
            self.gene = None
        else:
            self.gene = self.info["HGNC"]
    
    def passes_filters(self, filters):
        """Checks whether a VCF record passes user defined criteria.
        
        Args:
            filters: dict of filters
            
        Returns:
            boolean value for whether the record passes the filters
        """
        
        self.set_genotype()
        ddg2p = filters["HGNC"][1]
        
        track_variant = False
        if self.get_position() == "71295468":
            track_variant = True
        
        passes = True
        if self.fails_cnsolidate():
            passes = False
            if track_variant:
                print("failed CNSOLIDATE", self)
        elif self.fails_mad_ratio():
            passes = False
            if track_variant:
                print("failed mad ratio", self.info["MEANLR2"], self.info["MADL2R"])
        elif self.fails_wscore():
            passes = False
            if track_variant:
                print("failed wscore", self.info["WSCORE"])
        elif self.fails_callp():
            passes = False
            if track_variant:
                print("failed callp", self.info["CALLP"])
        elif self.fails_commmon_forwards():
            passes = False
            if track_variant:
                print("failed commonforwards", self.info["COMMONFORWARDS"])
        elif self.fails_meanlr2():
            passes = False
            if track_variant:
                print("failed meanlr2", self.info["MEANLR2"])
    
        return passes
    
    def fails_cnsolidate(self):
        """ checks if the CNV is passed by CNSOLIDATE
        """
        
        return "CNSOLIDATE" not in self.info
    
    def fails_mad_ratio(self):
        """ checks if the MAD ratio is too low
        """
        
        try:
            return abs(float(self.info["MEANLR2"])/float(self.info["MADL2R"])) < 15
        except ZeroDivisionError:
            return True
        
    def fails_wscore(self):
        """ checks if the WSCORE value is too low
        """
        
        return float(self.info["WSCORE"]) < 0.4
    
    def fails_callp(self):
        """ checks if the CALLP value is too high
        """
        
        return float(self.info["CALLP"]) > 0.01
    
    def fails_commmon_forwards(self):
        """ checks if the COMMONFORWARDS value is too high
        """
        
        return float(self.info["COMMONFORWARDS"]) > 0.8
    
    def fails_meanlr2(self):
        """ checks if the MEANLR2 value is out of bounds
        """
        
        if self.genotype == "DUP":
            return float(self.info["MEANLR2"]) < 0.36
        elif self.genotype == "DEL":
            return float(self.info["MEANLR2"]) > -0.41
        
        return False
        
    def is_het(self):
        # return False
        return self.genotype in self.alt_genotypes
    
    def is_hom_alt(self):
        # return False
        return self.genotype in self.alt_genotypes
    
    def is_hom_ref(self):
        # return True
        return self.genotype in self.ref_genotypes
    
    def is_not_ref(self):
        # return False
        return self.genotype in self.alt_genotypes
    
    def is_not_alt(self):
        # return True
        return self.genotype in self.ref_genotypes
