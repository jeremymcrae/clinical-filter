""" class for holding copy number information for a single individual
"""

from clinicalfilter.vcf_info import VcfInfo
from clinicalfilter.variant import Variant
from clinicalfilter.variant_cnv_acgh_filter import ACGH_CNV
from clinicalfilter.variant_cnv_exome_filter import ExomeCNV
from clinicalfilter.variant_cnv_breakdancer_filter import BreakdancerCNV

class CNV(Variant, VcfInfo):
    """  class to take CNV data for an individual, and
    """
    
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
            self.position = cnv_end
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
        
        self.ref_genotypes = set(["REF"])
        self.alt_genotypes = set(["DEL", "DUP"])
        
        if self.get_inheritance_type() == "YChrFemale":
            raise ValueError("cannot have CNV on female Y chromosome")
        
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
    
    def fix_gene_IDs(self, known_genes):
        """ find the genes that the CNV overlaps from a dict of known genes
        
        Sometimes the gene annotation for a CNV is incorrect - VEP annotated
        that the CNV overlaps a gene when other tools show there is not overlap.
        We correct for these by checking against a set of known genes 
        (currently the DDG2P set).
        
        Args:
            known_genes: dictionary of known genes, indexed by gene name, 
                or None
        """
        
        self.set_range()
        
        genes = []
        for gene in self.get_genes():
            # if the gene isn't in the DDG2P set we just include it as is, in
            # order to allow for non DDG2P variant analyses.
            # TODO: ideally we would match against all gencode positions
            if known_genes is None or gene not in known_genes:
                genes.append(gene)
            else:
                start = int(known_genes[gene]["start"])
                end = int(known_genes[gene]["end"])
                
                # only add the known gene if the DDG2P GENCODE positions
                # indicate that it overlaps with the CNV, otherwise exclude it.
                if int(self.start_position) <= end and int(self.end_position) >= start:
                    genes.append(gene)
        
        if len(genes) > 0:
            self.gene = ",".join(genes)
    
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
    
    def get_genes(self):
        """ split a gene string into list of gene names
        
        Returns:
            list of gene IDs
        """
        
        if self.gene == None:
            genes = []
        elif "," in self.gene:
            genes = self.gene.split(",")
        else:
            genes = [self.gene]
        
        return genes
    
    def passes_filters(self, filters):
        """Checks whether a VCF variant passes user defined criteria.
        
        Args:
            filters: filtering criteria (used for SNV filtering)
        
        Returns:
            boolean value for whether the variant passes the filters
        """
        
        # some CNVs are on female Y chrom, which give errors, fail those CNVs
        try:
            self.set_genotype()
        except ValueError:
            return False
        
        track_variant = False
        if self.get_position() == "186746704":
            track_variant = True
        
        passes = True
        # if self.fails_cnsolidate():
        #     passes = False
        #     if track_variant:
        #         print("failed CNSOLIDATE", self)
        # elif self.fails_mad_ratio():
        #     passes = False
        #     if track_variant:
        #         print("failed mad ratio", self.info["MEANLR2"], self.info["MADL2R"])
        # elif self.fails_wscore():
        #     passes = False
        #     if track_variant:
        #         print("failed wscore", self.info["WSCORE"])
        # elif self.fails_callp():
        #     passes = False
        #     if track_variant:
        #         print("failed callp", self.info["CALLP"])
        # elif self.fails_commmon_forwards():
        #     passes = False
        #     if track_variant:
        #         print("failed commonforwards", self.info["COMMONFORWARDS"])
        # elif self.fails_meanlr2():
        #     passes = False
        #     if track_variant:
        #         print("failed meanlr2", self.info["MEANLR2"])
        # elif self.fails_no_exons():
        #     passes = False
        #     if track_variant:
        #         print("failed no exons", self.info["NUMBEREXONS"])
        
        if "CONVEX" in self.info and "CNSOLIDATE" not in self.info:
            filt = ExomeCNV(self)
            passes = filt.filter_cnv(track_variant)
        elif "CNSOLIDATE" in self.info:
            filt = ACGH_CNV(self)
            passes = filt.filter_cnv(track_variant)
        elif "BREAKDANCER" in self.info:
            filt = BreakdancerCNV(self)
            passes = filt.filter_cnv(track_variant)
        else:
            if track_variant:
                 print("CNV is not an aCGH or exome CNV")
            passes = False
        
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
    
    def fails_no_exons(self):
        """ checks that the CNV overlaps at least one exon
        """
        
        return float(self.info["NUMBEREXONS"]) < 1
        
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
