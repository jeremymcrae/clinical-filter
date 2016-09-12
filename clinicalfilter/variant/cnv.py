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
from clinicalfilter.variant.cnv_acgh_filter import ACGH_CNV
from clinicalfilter.variant.cnv_exome_filter import ExomeCNV

class CNV(Variant):
    """  class for holding copy number information for a single individual
    """
    
    ref_genotypes = set(["REF"])
    alt_genotypes = set(["DEL", "DUP"])
    
    debug_chrom = None
    debug_pos = None
    
    @classmethod
    def set_debug(cls_obj, chrom, pos):
        cls_obj.debug_chrom = chrom
        cls_obj.debug_pos = pos
    
    def is_cnv(self):
        """ checks whether the variant is for a CNV
        """
        
        return True
    
    def set_genotype(self):
        """ sets the genotype of the variant
        """
        
        # ensure the inheritance type ("autosomal", "XChrMale" etc) is correct
        # for CNVs, since they can overlap allosomal and pseudoautosomal regions.
        self.set_inheritance_type(self.get_position(), self.is_male())
        start_inh = self.get_inheritance_type()
        
        self.set_inheritance_type(int(self.info["END"]), self.is_male())
        end_inh = self.get_inheritance_type()
        
        # CNVs that overlap allosomal and pseudoautosomal regions will have
        # different inheritance types for the range ends. Swap to allosomal.
        if start_inh != end_inh:
            # currently we are using the end inh type, so we only need to
            # swap to the start type if that is the allosomal end
            if start_inh != "autosomal":
                self.set_inheritance_type(start, self.is_male())
        
        if "CALLSOURCE" in self.info and self.info["CALLSOURCE"] == "EXOME":
            self.add_cns_state()
        
        if self.get_inheritance_type() == "YChrFemale" and '<REF>' not in self.alt_alleles:
            raise ValueError("cannot have CNV on female Y chromosome")
        
        if "<DUP>" in self.alt_alleles:
            self.genotype = "DUP"
        elif "<DEL>" in self.alt_alleles:
            self.genotype = "DEL"
        elif "<REF>" in self.alt_alleles:
            self.genotype = "REF"
        else:
            raise ValueError("unknown CNV allele code")
    
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
