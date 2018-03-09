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

import unittest
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.variant.cnv_exome_filter import ExomeCNV

from tests.utils import create_cnv

class TestExomeCnvPy(unittest.TestCase):
    """ test the Exome CNV filters
    """
    
    def setUp(self):
        """ define a default VcfInfo object
        """
        
        extra = "HGNC=OR5A1;CQ=missense_variant;CONVEX;RC50INTERNALFREQ=0.005;COMMONFORWARDS=0.000;MEANLR2=0.5;MADL2R=0.02"
        cnv = create_cnv('F', 'deNovo', extra_info=extra)
        
        self.var = ExomeCNV(cnv)
        
    
    def test_add_cns_state(self):
        """ test that add_cns_state() works correctly
        """
        
        # meanl2r > 0 gives cns of 3
        self.var.cnv.info["MEANLR2"] = "1"
        self.var.cnv.add_cns_state()
        self.assertEqual(self.var.cnv.info["CNS"], "3")
        
        self.var.cnv.info["MEANLR2"] = "0"
        self.var.cnv.add_cns_state()
        self.assertEqual(self.var.cnv.info["CNS"], "3")
        
        self.var.cnv.info["MEANLR2"] = "-0.1"
        self.var.cnv.add_cns_state()
        self.assertEqual(self.var.cnv.info["CNS"], "1")
        
        self.var.cnv.info["MEANLR2"] = "-2"
        self.var.cnv.add_cns_state()
        self.assertEqual(self.var.cnv.info["CNS"], "1")
        
        self.var.cnv.info["MEANLR2"] = "-2.1"
        self.var.cnv.add_cns_state()
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
        self.var.cnv.info["MADL2R"] = "0.01"
        self.assertFalse(self.var.fails_mad_ratio())
        
        # check that var fails when MAD ratio < 5
        self.var.cnv.info["MEANLR2"] = "0.09"
        self.var.cnv.info["MADL2R"] = "0.01"
        self.assertTrue(self.var.fails_mad_ratio())
        
        # check that var fails when trying to divide by zero
        self.var.cnv.info["MEANLR2"] = "0.2"
        self.var.cnv.info["MADL2R"] = "0"
        self.assertTrue(self.var.fails_mad_ratio())
    
    def test_fails_meanlr2(self):
        """ test that fails_meanlr2() works correctly
        """
        
        # below the DUP threshold
        self.var.cnv.genotype = 'DUP'
        self.var.cnv.info["MEANLR2"] = "0.399"
        self.assertTrue(self.var.fails_meanlr2())
        
        # at the DUP threshold
        self.var.cnv.info["MEANLR2"] = "0.4"
        self.assertFalse(self.var.fails_meanlr2())
        
        # above the DUP threshold
        self.var.cnv.info["MEANLR2"] = "0.41"
        self.assertFalse(self.var.fails_meanlr2())
        
        # above the DEL threshold
        self.var.cnv.genotype = 'DEL'
        self.var.cnv.info["MEANLR2"] = "-0.4999"
        self.assertTrue(self.var.fails_meanlr2())
        
        # at the DEL threshold
        self.var.cnv.info["MEANLR2"] = "-0.5"
        self.assertFalse(self.var.fails_meanlr2())
        
        # below the DEL threshold
        self.var.cnv.info["MEANLR2"] = "-0.51"
        self.assertFalse(self.var.fails_meanlr2())
    
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

    def test_fails_x_lr2(self):
        """test that fails_x_lr2() works correctly
        """
        self.var.cnv.chrom = "X"
        self.var.cnv.sum_x_lr2 = "-6000"
        self.assertTrue(self.var.fails_x_lr2())

        self.var.cnv.sum_x_lr2 = "8000"
        self.assertTrue(self.var.fails_x_lr2())

        self.var.cnv.sum_x_lr2 = "2000"
        self.assertFalse(self.var.fails_x_lr2())

    def test_fails_additional_filters(self):
        """test that fails_additional_filters() works correctly
        """
        self.var.cnv.genotype = "DUP"
        self.var.cnv.format["CIFER_INHERITANCE"] = "maternal_inh"
        self.var.cnv.info["MEANLR2"] = "-2"
        self.var.cnv.info["CONVEXSCORE"] = "10"
        self.var.cnv.info["MADL2R"] = "0.2"
        self.assertFalse(self.var.fails_additional_filters())

        self.var.cnv.genotype = "DEL"
        self.assertFalse(self.var.fails_additional_filters())

        self.var.cnv.format["CIFER_INHERITANCE"] = "not_inherited"
        self.assertTrue(self.var.fails_additional_filters())

        self.var.cnv.info["MADL2R"] = "0.1"
        self.assertTrue(self.var.fails_additional_filters())

        self.var.cnv.info["MEANLR2"] = "-1"
        self.assertFalse(self.var.fails_additional_filters())

        self.var.cnv.info["CONVEXSCORE"] = "20"
        self.assertFalse(self.var.fails_additional_filters())

        

if __name__ == '__main__':
    unittest.main()
