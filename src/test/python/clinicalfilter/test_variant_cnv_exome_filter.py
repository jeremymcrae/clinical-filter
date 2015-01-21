""" unit testing of the Exome CNV filtering class
"""

import unittest
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.variant.cnv_exome_filter import ExomeCNV


class TestExomeCnvPy(unittest.TestCase):
    """ test the Exome CNV filters
    """
    
    def setUp(self):
        """ define a default VcfInfo object
        """
        
        chrom = "1"
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "<DUP>"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        cnv = CNV(chrom, pos, snp_id, ref, alt, filt)
        self.var = ExomeCNV(cnv)
        
        info = "HGNC=TEST;HGNC_ALL=TEST,OR5A1;CQ=missense_variant;CONVEX;RC50INTERNALFREQ=0.005;COMMONFORWARDS=0.000;MEANLR2=0.5;MADL2R=0.02"
        format_keys = "inheritance:DP"
        sample_values = "deNovo:50"
        
        self.var.cnv.add_info(info)
        self.var.cnv.add_format(format_keys, sample_values)
        self.var.cnv.set_gender("F")
    
    def test_add_cns_state(self):
        """ test that add_cns_state() works correctly
        """
        
        # meanl2r > 0 gives cns of 3
        self.var.cnv.info["MEANLR2"] = "1"
        self.var.add_cns_state()
        self.assertEqual(self.var.cnv.info["CNS"], "3")
        
        self.var.cnv.info["MEANLR2"] = "0"
        self.var.add_cns_state()
        self.assertEqual(self.var.cnv.info["CNS"], "3")
        
        self.var.cnv.info["MEANLR2"] = "-0.1"
        self.var.add_cns_state()
        self.assertEqual(self.var.cnv.info["CNS"], "1")
        
        self.var.cnv.info["MEANLR2"] = "-2"
        self.var.add_cns_state()
        self.assertEqual(self.var.cnv.info["CNS"], "1")
        
        self.var.cnv.info["MEANLR2"] = "-2.1"
        self.var.add_cns_state()
        self.assertEqual(self.var.cnv.info["CNS"], "0")
    
    def test_fails_mad_ratio(self):
        """ test that fails_mad_ratio() works correctly
        """
        
        # check that var passes when MAD ratio > 5
        self.var.cnv.info["MEANLR2"] = "0.5"
        self.var.cnv.info["MADL2R"] = "0.02"
        self.assertFalse(self.var.fails_mad_ratio())
        
        # check that var passes when MAD ratio == 5
        self.var.cnv.info["MEANLR2"] = "0.1"
        self.var.cnv.info["MADL2R"] = "0.02"
        self.assertFalse(self.var.fails_mad_ratio())
        
        # check that var fails when MAD ratio < 5
        self.var.cnv.info["MEANLR2"] = "0.09"
        self.var.cnv.info["MADL2R"] = "0.02"
        self.assertTrue(self.var.fails_mad_ratio())
        
        # check that var fails when trying to divide by zero
        self.var.cnv.info["MEANLR2"] = "0.2"
        self.var.cnv.info["MADL2R"] = "0"
        self.assertTrue(self.var.fails_mad_ratio())
    
    def test_fails_population_frequency(self):
        """ test that fails_population_frequency() works correctly
        """
        
        # check that var passes when RC50INTERNALFREQ < 0.01
        self.var.cnv.info["RC50INTERNALFREQ"] = "0.005"
        self.assertFalse(self.var.fails_population_frequency())
        
        # check that var passes when RC50INTERNALFREQ == 0.1
        self.var.cnv.info["RC50INTERNALFREQ"] = "0.01"
        self.assertFalse(self.var.fails_population_frequency())
        
        # check that var fails when RC50INTERNALFREQ > 0.1
        self.var.cnv.info["RC50INTERNALFREQ"] = "0.02"
        self.assertTrue(self.var.fails_population_frequency())
    
    def test_fails_convex_score(self):
        """ test that fails_convex_score() works correctly
        """
        
        # check that var passes when CONVEXSCORE > 7
        self.var.cnv.info["CONVEXSCORE"] = "8"
        self.assertFalse(self.var.fails_convex_score())
        
        # check that var fails when CONVEXSCORE == 7
        self.var.cnv.info["CONVEXSCORE"] = "7"
        self.assertTrue(self.var.fails_convex_score())
        
        # check that var fails when CONVEXSCORE < 7
        self.var.cnv.info["CONVEXSCORE"] = "6"
        self.assertTrue(self.var.fails_convex_score())
    
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
    
    def test_fails_no_exons(self):
        """ test that fails_no_exons() works correctly
        """
        
        self.var.cnv.info["NUMBEREXONS"] = "1"
        self.assertFalse(self.var.fails_no_exons())

        self.var.cnv.info["NUMBEREXONS"] = "0"
        self.assertTrue(self.var.fails_no_exons())


if __name__ == '__main__':
    unittest.main()
