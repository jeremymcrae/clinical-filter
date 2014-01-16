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
        
        if self.alt_allele == "<DUP>":
            self.genotype = "DUP"
        elif self.alt_allele == "<DEL>":
            self.genotype = "DEL"
        else:
            self.alt_allele = "REF"
        
        self.set_range()
        self.set_reference_genotypes()
    
    def set_default_genotype(self):
        """ set a default genotype for individuals without one
        """
        
        self.genotype = "REF"
        self.set_range()
        self.set_reference_genotypes()
    
    def set_reference_genotypes(self):
        """ sets reference genotypes for homozygotes and heterozygotes, for checking against
        """
        
        self.ref_genotypes = set(["REF"])
        self.alt_genotypes = set(["DEL", "DUP"])
    
    def set_range(self):
        """ sets the range for the CNV
        """
        
        self.start_position = self.position
        
        if hasattr(self, "info"):
            self.end_position = self.info["END"]
        else:
            self.end_position = str(int(self.start_position) + 10000)
        
        self.range = (self.start_position, self.end_position)
        self.size = int(self.end_position) - int(self.start_position)
    
    def get_key(self):
        """ return a tuple to identify the variant
        """
        
        try:
            return (self.chrom, self.start_position, self.end_position)
        except AttributeError:
            self.set_range()
        
        return (self.chrom, self.start_position, self.end_position)
    
    def set_parental_statuses(self, mother_affected, father_affected):
        """ adds parental affected statuses for CNV filtering
        """
        
        self.father_affected = father_affected
        self.mother_affected = mother_affected
    
    def passes_filters(self, filters):
        """Checks whether a VCF record passes user defined criteria.
        
        Args:
            filters: dict of filters
            
        Returns:
            boolean value for whether the record passes the filters
        """
        
        return False
        
        self.set_genotype()
        
        track_variant = False
        if self.get_position() == "51135951":
            track_variant = True
        
        # select only acgh CNV - INFO.CNSOLIDATE [flag]
        if "CNSOLIDATE" not in self.info:
            if track_variant:
                print("failed CNSOLIDATE")
            return False
        
        try:
            if abs(float(self.info["MEANLR2"])/float(self.info["MADL2R"])) < 15:
                if track_variant:
                    print("failed MEANLR2/MADL2R, meanlr2:", self.info["MEANLR2"], "madl2r:", self.info["MADL2R"])
                return False
        except ZeroDivisionError:
            if track_variant:
                print("had zero division")
            return False
        
        if float(self.info["WSCORE"]) < 0.4:
            if track_variant:
                print("failed WSCORE", self.info["WSCORE"])
            return False
        
        if float(self.info["CALLP"]) > 0.01:
            if track_variant:
                print("failed CALLP", self.info["CALLP"])
            return False
        
        if float(self.info["COMMONFORWARDS"]) > 0.8:
            if track_variant:
                print("failed COMMONFORWARDS", self.info["COMMONFORWARDS"])
            return False
       
        if self.genotype == "DUP":
            if float(self.info["MEANLR2"]) < 0.36:
                if track_variant:
                    print("failed MEANLR2 for DUP", self.info["MEANLR2"])
                return False
        elif self.genotype == "DEL":
            if float(self.info["MEANLR2"]) > -0.41:
                if track_variant:
                    print("failed MEANLR2 for DEL", self.info["MEANLR2"])
                return False
        
        if self.gene is None or self.gene not in filters["HGNC"][1]:
            inh = self.format["INHERITANCE"]
            
            unknown_inh = set(["uknown", "inconclusive", "NA", "noCNVfound", \
                "inconclusiveDuo"])
            
            if inh == "deNovo" \
               or (inh == "paternal" and self.father_affected) \
               or (inh == "maternal" and self.mother_affected) \
               or (inh == "biparental" and \
                (self.father_affected or self.mother_affected)):
                if self.genotype == "DEL":
                    if float(self.info["SVLEN"]) < 100000:
                        if track_variant:
                            print("failed SVLEN for DEL", self.info["SVLEN"])
                        return False
                if self.genotype == "DUP":
                    if float(self.info["SVLEN"]) < 250000:
                        if track_variant:
                            print("failed SVLEN for DUP", self.info["SVLEN"])
                        return False
            
            elif inh in unknown_inh:
                if float(self.info["SVLEN"]) < 500000:
                    if track_variant:
                        print("failed SVLEN for unknown", self.info["SVLEN"])
                    return False
            else:
                if track_variant:
                    print("failed non-DDG2P inheritance", self, inh)
                return False
        
        elif self.gene in filters["HGNC"][1]:
            dd_gene_type = filters["HGNC"][1][self.gene]["confirmed_status"]
            
            if "Both DD and IF" in dd_gene_type:
                pass
            elif "Confirmed DD Gene" in dd_gene_type or "Probable DD gene" in dd_gene_type:
                required_chrom, required_copy_number, required_mechanisms, inheritance = self.set_CNV_filter_values(filters)
                
                if required_chrom == "X" and self.get_chrom() != "X":
                    if track_variant:
                        print("failed required chrom")
                    return False
                
                if self.info["CNS"] not in required_copy_number:
                    if track_variant:
                        print(filters["HGNC"][1][self.gene])
                        print(required_chrom, required_copy_number, required_mechanisms, inheritance)
                        print("failed required copy number state", self, self.info["CNS"], "not in", required_copy_number)
                    return False
                
                if len(filters["HGNC"][1][self.gene]["inheritance"][inheritance] & required_mechanisms) == 0:
                    if track_variant:
                        print("failed required inheritance mechanism", filters["HGNC"][1][self.gene]["inheritance"][inheritance], "not in", required_mechanisms)
                    return False
            else:
                print(self, "failed some DDG2P CNV criteria")
        
        if self.gene in filters["HGNC"][1]:
            print("filtered", self, filters["HGNC"][1][self.gene])
        else:
            print("filtered", self)

        return True
    
    def set_CNV_filter_values(self, filters):
        """ create CNV filter values dependent on DDG2P inheritance mode
        """
        
        dd_gene_mode = filters["HGNC"][1][self.gene]["inheritance"]
        
        if "Biallelic" in dd_gene_mode:
            required_chrom = "all"
            required_copy_number = set(["0"])
            required_mechanisms = set(["Uncertain", "Loss of function", "Dominant negative"])
            inheritance = "Biallelic"
        elif "Monoallelic" in dd_gene_mode:
            required_chrom = "all"
            required_copy_number = (["0", "1", "3"])
            required_mechanisms = set(["Uncertain", "Loss of function", "Dominant negative", "Increased gene dosage"])
            inheritance = "Monoallelic"
        elif "X-linked dominant" in dd_gene_mode:
            required_chrom = "X"
            required_copy_number = set(["0", "1", "3"])
            required_mechanisms = set(["Uncertain", "Loss of function", "Dominant negative", "Increased gene dosage"])
            inheritance = "X-linked dominant"
        elif "Hemizygous" in dd_gene_mode and self.is_male():
            required_chrom = "X"
            required_copy_number = set(["0", "1", "3"])
            required_mechanisms = set(["Uncertain", "Loss of function", "Dominant negative", "Increased gene dosage"])
            inheritance = "Hemizygous"
        elif "Hemizygous" in dd_gene_mode and self.is_female():
            required_chrom = "X"
            required_copy_number = set(["3"])
            required_mechanisms = set(["Increased gene dosage"])
            inheritance = "Hemizygous"
        else:
            # other inheritance modes of "Mosaic", or "Digenic" can be ignored
            # by using impossible criteria
            required_chrom = "GGG"
            required_copy_number = set(["999"])
            required_mechanisms = set(["Nothing"])
            inheritance = list(set(dd_gene_mode))[0]
            
        return required_chrom, required_copy_number, required_mechanisms, inheritance
    
    def is_het(self):
        # return False
        if self.genotype in self.alt_genotypes:
            return True
        else:
            return False
    
    def is_hom_alt(self):
        # return False
        if self.genotype in self.alt_genotypes:
            return True
        else:
            return False
    
    def is_hom_ref(self):
        # return True
        if self.genotype in self.ref_genotypes:
            return True
        else:
            return False
    
    def is_not_ref(self):
        # return False
        if self.genotype in self.alt_genotypes:
            return True
        else:
            return False
    
    def is_not_alt(self):
        # return True
        if self.genotype in self.ref_genotypes:
            return True
        else:
            return False
