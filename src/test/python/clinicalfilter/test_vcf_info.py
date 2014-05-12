""" unit testing of the VcfInfo class
"""

import unittest
from clinicalfilter.variant_snv import SNV

class TestVcfInfoPy(unittest.TestCase):
    """
    """
    
    def setUp(self):
        """ define a default VcfInfo object
        """
        
        chrom = "1"
        pos = "15000000"
        snp_id = "CM00001"
        ref = "A"
        alt = "G"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        self.var = SNV(chrom, pos, snp_id, ref, alt, filt)
        
        self.default_info = "HGNC=ATRX;CQ=missense_variant;random_tag"
        self.default_tags = {"gene": ["HGNC", "VGN", "GN"], "consequence": \
            ["VCQ", "CQ"], "MAX_MAF": ["AF_MAX", "MAX_AF", "1000G_AF", \
            "1KG_AF", "AF_AFR", "AF_AMR", "AF_ASN", "AF_EUR", \
            "UK10K_cohort_AF", "ESP_AF", "DDD_AF", "ASN_AF", "AFR_AF", \
            "EUR_AF", "AMR_AF"], "transcript": ["ENSEMBL_TRANSCRIPT", "ENST"]}
        
        # here are the default filtering criteria, as loaded into python
        self.default_filters = {"FILTER": ("list", set(["PASS", "."])), \
            "VCQ": ("list", set(["ESSENTIAL_SPLICE_SITE", "STOP_GAINED", \
                "COMPLEX_INDEL", "FRAMESHIFT_CODING", "NON_SYNONYMOUS_CODING", \
                "STOP_LOST", "transcript_ablation", "splice_donor_variant", \
                "splice_acceptor_variant", "frameshift_variant", \
                "initiator_codon_variant", "inframe_insertion", \
                "inframe_deletion", "missense_variant", \
                "transcript_amplification", "stop_gained", "stop_lost", \
                "coding_sequence_variant"])), \
            "MAX_MAF": ("smaller_than", 0.01)}
        
        # add all the populations that have minor allele frequency data in the 
        # VCF files into the filtering criteria
        for pop in self.default_tags["MAX_MAF"]:
            self.default_filters[pop] = self.default_filters["MAX_MAF"]
        
        # add all the possible consequence tags to the filter
        for cq in self.default_tags["consequence"]:
            self.default_filters[cq] = self.default_filters["VCQ"]
        
        self.var.add_info(self.default_info, self.default_tags)
    
    def test_add_gene_from_info(self):
        """ test that test_add_gene_from_info() works correctly
        """
        
        # check for when a HGNC key exists
        self.var.info["HGNC"] = "A"
        self.var.add_gene_from_info()
        self.assertEqual(self.var.gene, "A")
        
        # check for when a HGNC key doesn't exist
        del self.var.info["HGNC"]
        self.var.add_gene_from_info()
        self.assertIsNone(self.var.gene)
    
    def test_is_lof(self):
        """ test that is_lof() works correctly
        """
        
        # check that known LOF consensequence return True
        self.var.info["CQ"] = "stop_gained"
        self.assertTrue(self.var.is_lof())
        
        # check that known non-LOF consensequence returns False
        self.var.info["CQ"] = "missense_variant"
        self.assertFalse(self.var.is_lof())
        
        # check that null values return False
        self.var.info["CQ"] = None
        self.assertFalse(self.var.is_lof())
    
    def test_get_number(self):
        """ tests that number conversion works as expected
        """
        
        # single number returns that number
        self.assertEqual(self.var.get_number("1"), 1)
        
        # two numbers return one number
        self.assertEqual(self.var.get_number("1,1"), 1)
        
        # two numbers return the first number
        self.assertEqual(self.var.get_number("1,2"), 1)
        
        # number and string return the number
        self.assertEqual(self.var.get_number("a,1"), 1)
        
        # single string value returns that string
        self.assertEqual(self.var.get_number("a"), "a")
        
        # multiple string values return the first string value
        self.assertEqual(self.var.get_number("a,b"), "b")
    
    def test_is_number(self):
        """ tests that we can check if a value represents a number
        """
        
        self.assertEqual(self.var.is_number(None), False)
        self.assertEqual(self.var.is_number("5"), True)
        self.assertEqual(self.var.is_number("a"), False)
    
    def test_find_max_allele_frequency(self):
        """ test if the MAF finder operates correctly
        """
        
        pops = self.default_tags["MAX_MAF"]
        
        # check for var without recorded MAF
        self.assertEqual(self.var.find_max_allele_frequency(pops), "NA")
        
        # check for single population
        self.var.info["AF_MAX"] = "0.005"
        self.assertEqual(self.var.find_max_allele_frequency(pops), "0.005")
        
        # check for two populations
        self.var.info["AF_AFR"] = "0.01"
        self.assertEqual(self.var.find_max_allele_frequency(pops), "0.01")
        
        # check for all populations
        for pop in pops:
            self.var.info["AF_AFR"] = "0.05"
        self.assertEqual(self.var.find_max_allele_frequency(pops), "0.05")
    
    def test_passes_default_filters(self):
        """ test that different variants pass or fail the VcfInfo filters
        """
        
        # check that a default variant passes the filters
        self.assertTrue(self.var.passes_filters(self.default_filters))
    
    def test_passes_alternate_filter_string(self):
        """ test that the alternate permitted FILTER string also passes
        """
        
        # check that the alternate FILTER value passes
        self.var.info["FILTER"] = "."
        self.assertTrue(self.var.passes_filters(self.default_filters))
        
    def test_passes_filters_low_maf(self):
        """ tests that low MAF values pass the filters
        """
        
        # check that low MAF values pass the filters
        pops = self.default_tags["MAX_MAF"]
        for pop in pops:
            self.var.info[pop] = "0.005"
            self.assertTrue(self.var.passes_filters(self.default_filters))
            
            # and check that MAF on the threshold still pass
            self.var.info[pop] = "0.01"
            self.assertTrue(self.var.passes_filters(self.default_filters))
    
    def test_out_of_range_maf(self):
        """ check that MAF outside 0-1 still pass or fail correctly
        """
        
        self.var.info["AF_AFR"] = "-1"
        self.assertTrue(self.var.passes_filters(self.default_filters))
      
        self.var.info["AF_AFR"] = "100"
        self.assertFalse(self.var.passes_filters(self.default_filters))
      
    def test_fails_filters_high_maf(self):
        """ test that variants with high MAF fail the filtering
        """
        
        pops = self.default_tags["MAX_MAF"]
        
        # check th
        for pop in pops:
            var = self.var
            var.info[pop] = "0.0101"
            self.assertFalse(var.passes_filters(self.default_filters))
    
    def test_passes_consequence_filter(self):
        """ check all the consequence values that should pass
        """
        
        vep_passing = ["transcript_ablation", "splice_donor_variant", \
            "splice_acceptor_variant", "frameshift_variant", \
            "initiator_codon_variant", "inframe_insertion", "inframe_deletion",\
            "missense_variant", "transcript_amplification", "stop_gained",\
            "stop_lost", "coding_sequence_variant"]
        
        # check all the passing consequences
        for cq in vep_passing:
            self.var.info["CQ"] = cq
            self.assertTrue(self.var.passes_filters(self.default_filters))
            
    def test_fails_consequence_filter(self):
        """ check all the consequence values that should fail
        """
        
        vep_failing = ["splice_region_variant", \
            "incomplete_terminal_codon_variant", "synonymous_variant", \
            "stop_retained_variant", "mature_miRNA_variant", \
            "5_prime_UTR_variant", "3_prime_UTR_variant", \
            "non_coding_exon_variant", "nc_transcript_variant", \
            "intron_variant", "NMD_transcript_variant", \
            "upstream_gene_variant", "downstream_gene_variant", \
            "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", \
            "regulatory_region_variant", "regulatory_region_ablation", \
            "regulatory_region_amplification", "feature_elongation", \
            "feature_truncation", "intergenic_variant"]
        
        # check all the failing consequences
        for cq in vep_failing:
            self.var.info["CQ"] = cq
            self.assertFalse(self.var.passes_filters(self.default_filters))
    
    def test_passes_list(self):
        """ tests that passes_list() operates correctly
        """
        #  check for value in list
        self.assertTrue(self.var.passes_list("1", ["1", "2", "3"]))
        
        # check for value not in list
        self.assertFalse(self.var.passes_list("4", ["1", "2", "3"]))
        
        # check for None value
        self.assertFalse(self.var.passes_list(None, ["1", "2", "3"]))
        
        # check for zero length list
        self.assertFalse(self.var.passes_list("1", []))
        
        # check for None list
        self.assertFalse(self.var.passes_list("1", None))
    
    def test_passes_smaller_than(self):
        """ tests that passes_smaller_than() operates correctly
        """
        
        # check for smaller string value
        self.assertTrue(self.var.passes_smaller_than("0", 1))
        
        # check for larger string value
        self.assertFalse(self.var.passes_smaller_than("2", 1))
        
        # check for nonconvertible string
        self.assertTrue(self.var.passes_smaller_than("a", 1))
        
        # check for None value
        self.assertTrue(self.var.passes_smaller_than(None, 1))


if __name__ == '__main__':
    unittest.main()


