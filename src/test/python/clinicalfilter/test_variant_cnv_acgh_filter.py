""" unit testing of the ACGH CNV filtering class
"""

import unittest
from clinicalfilter.variant_cnv import CNV
from clinicalfilter.variant_cnv_acgh_filter import ACGH_CNV


class TestAcghCnvPy(unittest.TestCase):
    """ test the ACGH CNV filters
    """
    
    def setUp(self):
        """ define a default VcfInfo object
        """
        
        chrom = "1"
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "<DUP>"
        qual = "50"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        cnv = CNV(chrom, pos, snp_id, ref, alt, qual, filt)
        self.var = ACGH_CNV(cnv)
        
        tags = {"gene": ["HGNC", "VGN", "GN"], "consequence": \
            ["VCQ", "CQ"]}
        
        info = "HGNC=TEST;HGNC_ALL=TEST,OR5A1;CQ=missense_variant;CNSOLIDATE;WSCORE=0.5;CALLP=0.000;COMMONFORWARDS=0.000;MEANLR2=0.5;MADL2R=0.02;END=16000000;SVLEN=1000000"
        format_keys = "inheritance:DP"
        sample_values = "deNovo:50"
        
        self.var.cnv.add_info(info, tags)
        self.var.cnv.add_format(format_keys, sample_values)
        self.var.cnv.set_gender("F")
    
    def test_fails_mad_ratio(self):
        """ test that fails_mad_ratio() works correctly
        """
        
        # check that var passes when MAD ratio > 15
        self.var.cnv.info["MEANLR2"] = "0.5"
        self.var.cnv.info["MADL2R"] = "0.02"
        self.assertFalse(self.var.fails_mad_ratio())
        
        # check that var passes when MAD ratio == 15
        self.var.cnv.info["MEANLR2"] = "0.3"
        self.var.cnv.info["MADL2R"] = "0.02"
        self.assertFalse(self.var.fails_mad_ratio())
        
        # check that var fails when MAD ratio < 15
        self.var.cnv.info["MEANLR2"] = "0.2"
        self.var.cnv.info["MADL2R"] = "0.02"
        self.assertTrue(self.var.fails_mad_ratio())
        
        # check that var fails when trying to divide by zero
        self.var.cnv.info["MEANLR2"] = "0.2"
        self.var.cnv.info["MADL2R"] = "0"
        self.assertTrue(self.var.fails_mad_ratio())
    
    def test_fails_wscore(self):
        """ test that fails_wscore() works correctly
        """
        
        # check that var passes when WSCORE > 0.4
        self.var.cnv.info["WSCORE"] = "0.5"
        self.assertFalse(self.var.fails_wscore())
        
        # check that var passes when WSCORE == 0.4
        self.var.cnv.info["WSCORE"] = "0.4"
        self.assertFalse(self.var.fails_wscore())
        
        # check that var fails when WSCORE > 0.4
        self.var.cnv.info["WSCORE"] = "0.399"
        self.assertTrue(self.var.fails_wscore())
    
    def test_fails_callp(self):
        """ test that fails_callp() works correctly
        """
        
        # check that var passes when CALLP < 0.01
        self.var.cnv.info["CALLP"] = "0.0"
        self.assertFalse(self.var.fails_callp())
        
        # check that var passes when CALLP == 0.01
        self.var.cnv.info["CALLP"] = "0.01"
        self.assertFalse(self.var.fails_callp())
        
        # check that var fails when CALLP > 0.01
        self.var.cnv.info["CALLP"] = "0.0101"
        self.assertTrue(self.var.fails_callp())
    
    def test_fails_commmon_forwards(self):
        """ test that fails_commmon_forwards() works correctly
        """
        
        # check that var passes when COMMONFORWARDS < 0.8
        self.var.cnv.info["COMMONFORWARDS"] = "0.0"
        self.assertFalse(self.var.fails_commmon_forwards())
        
        # check that var fails when COMMONFORWARDS == 0.8
        self.var.cnv.info["COMMONFORWARDS"] = "0.8"
        self.assertFalse(self.var.fails_commmon_forwards())
        
        # check that var fails when COMMONFORWARDS > 0.8
        self.var.cnv.info["COMMONFORWARDS"] = "0.9"
        self.assertTrue(self.var.fails_commmon_forwards())
    
    def test_fails_mean_lr2_dup(self):
        """ test that fails_mean_lr2() works correctly on duplications
        """
        
        # set the var as a duplication
        self.var.cnv.alt_allele = "<DUP>"
        self.var.cnv.set_genotype()
        
        # check that dup passes with MEANLR2 > 0.36
        self.var.cnv.info["MEANLR2"] = "0.4"
        self.assertFalse(self.var.fails_meanlr2())
        
        # check that dup passes with MEANLR2 == 0.36
        self.var.cnv.info["MEANLR2"] = "0.36"
        self.assertFalse(self.var.fails_meanlr2())
        
        # check that dup passes with MEANLR2 < 0.36
        self.var.cnv.info["MEANLR2"] = "0.359"
        self.assertTrue(self.var.fails_meanlr2())
    
    def test_fails_mean_lr2_del(self):
        """ test that fails_mean_lr2() works correctly on deletions
        """
        # set the var as a deletion
        self.var.cnv.alt_allele = "<DEL>"
        self.var.cnv.set_genotype()
        
        # check that del passes with MEANLR2 < -0.41
        self.var.cnv.info["MEANLR2"] = "-0.5"
        self.assertFalse(self.var.fails_meanlr2())
        
        # check that del passes with MEANLR2 == -0.41
        self.var.cnv.info["MEANLR2"] = "-0.41"
        self.assertFalse(self.var.fails_meanlr2())
        
        # check that del passes with MEANLR2 > -0.41
        self.var.cnv.info["MEANLR2"] = "-0.409"
        self.assertTrue(self.var.fails_meanlr2())
    
    def test_fails_no_exons(self):
        """ test that fails_no_exons() works correctly
        """
        
        self.var.cnv.info["NUMBEREXONS"] = "1"
        self.assertFalse(self.var.fails_no_exons())

        self.var.cnv.info["NUMBEREXONS"] = "0"
        self.assertTrue(self.var.fails_no_exons())


if __name__ == '__main__':
    unittest.main()


