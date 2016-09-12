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
from clinicalfilter.variant.cnv_acgh_filter import ACGH_CNV

from tests.utils import create_cnv

class TestAcghCnvPy(unittest.TestCase):
    """ test the ACGH CNV filters
    """
    
    def setUp(self):
        """ define a default VcfInfo object
        """
        
        extra = 'OR5A1;CNSOLIDATE;WSCORE=0.5;CALLP=0.000;COMMONFORWARDS=0.000;MEANLR2=0.5;MADL2R=0.02;'
        cnv = create_cnv('F', 'deNovo', extra_info=extra)
        
        self.var = ACGH_CNV(cnv)
    
    def test_fails_mad_ratio(self):
        """ test that fails_mad_ratio() works correctly.
        """
        
        # check that var passes when MAD ratio > 0
        self.var.cnv.info["MEANLR2"] = "0.5"
        self.var.cnv.info["MADL2R"] = "0.02"
        self.assertFalse(self.var.fails_mad_ratio())
        
        # check that var passes when MAD ratio == 0
        self.var.cnv.info["MEANLR2"] = "0.3"
        self.var.cnv.info["MADL2R"] = float("inf")
        self.assertFalse(self.var.fails_mad_ratio())
        
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
        self.var.cnv.alt_alleles = ["<DUP>"]
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
        self.var.cnv.alt_alleles = ["<DEL>"]
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
    
    def test_fails_frequency(self):
        """ test that fails_frequency() works correctly
        """
        
        # a low population frequency will pass
        self.var.cnv.info["ACGH_RC_FREQ50"] = "0.01"
        self.assertFalse(self.var.fails_frequency())
        
        # a high population frequency will fail
        self.var.cnv.info["ACGH_RC_FREQ50"] = "0.011"
        self.assertTrue(self.var.fails_frequency())
        
        # if population frequency information is unavailable, it should pass,
        # since this suggests the frequency is 0.
        del self.var.cnv.info["ACGH_RC_FREQ50"]
        self.assertFalse(self.var.fails_frequency())
    
    def test_fails_cifer_inh(self):
        """ test that fails_cifer_inh() works correctly
        """
        
        # CNVs annotated as not_inherited, or inherited will pass
        self.var.cnv.format["CIFER_INHERITANCE"] = "not_inherited"
        self.assertFalse(self.var.fails_cifer_inh())
        
        # CNVs annotated as false_positive will fail
        self.var.cnv.format["CIFER_INHERITANCE"] = "false_positive"
        self.assertTrue(self.var.fails_cifer_inh())
