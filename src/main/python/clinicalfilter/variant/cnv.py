""" class for holding copy number information for a single individual
"""

from clinicalfilter.variant.info import VariantInfo
from clinicalfilter.variant.variant import Variant
from clinicalfilter.variant.cnv_acgh_filter import ACGH_CNV
from clinicalfilter.variant.cnv_exome_filter import ExomeCNV

class CNV(Variant, VariantInfo):
    """  class to take CNV data for an individual, and
    """
    
    ref_genotypes = set(["REF"])
    alt_genotypes = set(["DEL", "DUP"])
    
    def is_cnv(self):
        """ checks whether the variant is for a CNV
        """
        
        return True
    
    def set_genotype(self):
        """ sets the genotype of the variant
        """
        
        # make sure the inheritance type ("autosomal", "XChrMale" etc) is
        # set correctly for allosomal CNVs, since they may lie across both
        # allosomal and pseudoautosomal regions.
        if self.get_chrom() in ["chrX", "ChrX", "X", "chrY", "ChrY", "Y"]:
            cnv_start = self.get_position()
            cnv_start_inh = self.get_inheritance_type()
            
            cnv_end = self.info["END"]
            self.position = int(cnv_end)
            self.set_inheritance_type()
            cnv_end_inh = self.get_inheritance_type()
            
            # restore the CNVs initial position
            self.position = cnv_start
            
            # if the start and end positions have different inheritance types,
            # swap ourselves over to the allosomal inheritance type
            if cnv_start_inh != cnv_end_inh:
                # currently we are using the end inh type, so we only need to
                # swap to the start type if that is the allosomal end
                if cnv_start_inh != "autosomal":
                    self.set_inheritance_type()
        
        if self.get_inheritance_type() == "YChrFemale":
            raise ValueError("cannot have CNV on female Y chromosome")
        
        if self.alt_allele == "<DUP>":
            self.genotype = "DUP"
        elif self.alt_allele == "<DEL>":
            self.genotype = "DEL"
        else:
            raise ValueError("unknown CNV allele code")
    
    def set_default_genotype(self):
        """ set a default genotype for individuals without one
        """
        
        self.genotype = "REF"
    
    def get_key(self):
        """ return a tuple to identify the variant
        """
        
        start, end = self.get_range()
        
        return (self.get_chrom(), start, end)
    
    def fix_gene_IDs(self):
        """ find the genes that the CNV overlaps from a dict of known genes
        
        Sometimes the gene annotation for a CNV is incorrect - VEP annotated
        that the CNV overlaps a gene when other tools show there is not overlap.
        We correct for these by checking against a set of known genes
        (currently the DDG2P set).
        """
        
        (start, end) = self.get_range()
        
        genes = []
        for gene in self.get_genes():
            # if the gene isn't in the DDG2P set we just include it as is, in
            # order to allow for non DDG2P variant analyses.
            # TODO: ideally we would match against all gencode positions
            if self.known_genes is None or gene not in self.known_genes:
                genes.append(gene)
            else:
                gene_start = self.known_genes[gene]["start"]
                gene_end = self.known_genes[gene]["end"]
                
                # only add the known gene if the DDG2P GENCODE positions
                # indicate that it overlaps with the CNV, otherwise exclude it.
                if start <= gene_end and end >= gene_start:
                    genes.append(gene)
                else:
                    genes.append(".")
        
        self.genes = genes
    
    def passes_filters(self):
        """Checks whether a VCF variant passes user defined criteria.
        
        Returns:
            boolean value for whether the variant passes the filters
        """
        
        # some CNVs are on female Y chrom, which give errors, fail those CNVs
        try:
            self.set_genotype()
        except ValueError:
            return False
        
        # we rely on the CALLSOURCE field to inform us what the CNV has been
        # called by. Raise an error if this is not present.
        assert "CALLSOURCE" in self.info
        
        track_variant = False
        if self.get_chrom() == self.debug_chrom and self.get_position() == self.debug_pos:
            track_variant = True
        
        passes = True
        if "aCGH" in self.info["CALLSOURCE"]:
            filt = ACGH_CNV(self)
            passes = filt.filter_cnv(track_variant)
        elif "EXOME" in self.info["CALLSOURCE"]:
            # currently return false for all exome-only CNVs, undergoing testing
            filt = ExomeCNV(self)
            self.add_cns_state()
            # passes = filt.filter_cnv(track_variant)
            passes = False
        else:
            if track_variant:
                print("CNV is not an aCGH or exome CNV")
            passes = False
        
        return passes
    
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
    
    def add_cns_state(self):
        """ determines the CNS value from MEANLR2 values
        """
        
        if float(self.info["MEANLR2"]) >= 0:
            self.info["CNS"] = "3"
        elif 0 > float(self.info["MEANLR2"]) >= -2:
            self.info["CNS"] = "1"
        elif -2 > float(self.info["MEANLR2"]):
            self.info["CNS"] = "0"
        else:
            raise ValueError("Shouldn't reach here")
