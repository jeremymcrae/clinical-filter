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
from clinicalfilter.variant.variant import Variant
from clinicalfilter.variant.info import Info


class TestVariantPy(unittest.TestCase):
    """ test Variant objects and their methods
    """
    
    def setUp(self):
        """ define a default variant (currently a SNV)
        """
        chrom = "1"
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "G"
        qual = "1000"
        filt = "PASS"
        info = 'AC=1;HGNC=TEST'
        keys = 'GT:DP:AD'
        sample = '0/1:40:10,10'
        
        self.var = Variant(chrom, pos, snp_id, ref, alt, qual, filt, info, keys, sample)
    
    def test_set_gender_unknown(self):
        """ tests set_gender(), and its implications on the inheritance type
        """
        
        # raise error for unknown gender
        with self.assertRaises(ValueError):
            self.var._set_gender("unknown")
    
    def test_set_gender_autosomal(self):
        """test gender on autosomal chroms
        """
        
        # check all the legitimate gender codes
        gender_codes = ["1", "m", "M", "male", "2", "f", "F", "female"]
        for gender in gender_codes:
            self.var._set_gender(gender)
            self.assertEqual(self.var.get_inheritance_type(), "autosomal")
    
    def test_set_gender_allosomal(self):
        """ test gender on the X chromosome
        """
        
        # check the allosomal chroms
        self.var.chrom = "X"
        self.var._set_gender("male")
        self.assertEqual(self.var.get_inheritance_type(), "XChrMale")
        
        self.var._set_gender("female")
        self.assertEqual(self.var.get_inheritance_type(), "XChrFemale")
        
         # check the allosomal chroms
        self.var.chrom = "Y"
        self.var._set_gender("male")
        self.assertEqual(self.var.get_inheritance_type(), "YChrMale")
        
        self.var._set_gender("female")
        self.assertEqual(self.var.get_inheritance_type(), "YChrFemale")
    
    def test_set_gender_pseudoautosomal(self):
        """ test gender on a pseudoautosomal region
        """
        
        # now check a variant in the pseudoautosomal regions
        self.var.chrom = "X"
        self.var.position = 2699510
        self.var._set_gender("male")
        self.assertEqual(self.var.get_inheritance_type(), "autosomal")
        
        self.var._set_gender("female")
        self.assertEqual(self.var.get_inheritance_type(), "autosomal")
    
    def test_get_chrom(self):
        """ tests that the chrom returns correctly
        """
        
        # check that the chrom is the default chrom
        self.assertEqual(self.var.get_chrom(), "1")
        
        # check that we can change the chrom and it still returns correctly
        self.var.chrom = "X"
        self.assertEqual(self.var.get_chrom(), "X")
    
    def test_get_position(self):
        """ tests that the position returns correctly
        """
        
        # check the default position
        self.assertEqual(self.var.get_position(), 15000000)
        
        # check that we can change the chrom and it still returns correctly
        self.var.position = 123456789
        self.assertEqual(self.var.get_position(), 123456789)
    
    def test_get_mutation_id(self):
        """ test that the mutation ID is parsed correctly
        """
        var = self.var
        
        # check a null ID value
        var.set_mutation_id(".")
        self.assertEqual(var.get_mutation_id(), "NA")
        
        # check a rsID value
        var.set_mutation_id("rs12546")
        self.assertEqual(var.get_mutation_id(), "NA")
        
        # check a rsID and a mutation ID
        var.set_mutation_id("rs12546&CM0001")
        self.assertEqual(var.get_mutation_id(), "CM0001")
        
        # check multiple mutation IDs
        var.set_mutation_id("CM0001&CM0002")
        self.assertEqual(var.get_mutation_id(), "CM0001,CM0002")
    
    def test_get_vcf_line(self):
        """ tests that the vcf line string returns correctly
        """
        
        vcf_line = ["1", "15000000", ".", "A", "G", "50", "PASS", "AB=0.41;AC=1;AN=2", "GT:gatk_PL:GQ", "0/1:736,0,356:99"]
        self.var.add_vcf_line(vcf_line)
        self.assertEqual(self.var.get_vcf_line(), vcf_line)
    
    def test_get_low_depth_alleles(self):
        ''' test that get_low_depth_alleles() works correctly
        '''
        
        # check with a single allele whre it is non-zero
        self.var.info = Info('AC=1')
        alts = ('C', )
        self.assertEqual(self.var.get_low_depth_alleles('G', alts), [])
        
        # check with a single allele with zero count
        self.var.info = Info('AC=0')
        alts = ('C', )
        self.assertEqual(self.var.get_low_depth_alleles('G', alts), ['C'])
        
        # check with multiallelic, where both are nonzero
        self.var.info = Info('AC=1,1')
        self.var.format = {'AD': '5,10,10'}
        alts = ('C', 'G')
        self.assertEqual(self.var.get_low_depth_alleles('G', alts), [])
        
        # check with multiallelic, where one a has zero count
        self.var.info = Info('AC=1,0')
        self.var.format = {'AD': '5,10,10'}
        alts = ('C', 'G')
        self.assertEqual(self.var.get_low_depth_alleles('G', alts), ['G'])
    
    def test_get_low_depth_alleles_bad_indel(self):
        ''' test that get_low_depth_alleles() works for indels with depth=1
        '''
        # check with multiallelic, where all should pass, since none are indels
        self.var.info = Info('AC=1,1')
        self.var.format = {'AD': '5,10,1'}
        alts = ('C', 'G')
        self.assertEqual(self.var.get_low_depth_alleles('G', alts), [])
        
        # but if we have an indel alt, and the corresponding depth is bad, fail
        alts = ('C', 'GG')
        self.assertEqual(self.var.get_low_depth_alleles('G', alts), ['GG'])
    
    # TODO: check add_format


if __name__ == '__main__':
    unittest.main()
