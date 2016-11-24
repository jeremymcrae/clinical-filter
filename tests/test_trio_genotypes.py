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
import copy

from clinicalfilter.ped import Family
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.variant.snv import SNV
from clinicalfilter.trio_genotypes import TrioGenotypes

from tests.utils import create_snv

class TestTrioGenotypesPy(unittest.TestCase):
    """ test the TrioGenotypes class
    """
    
    def setUp(self):
        """ define a family and variant, and start the Inheritance class
        """
        
        # generate a test family
        sex = "F"
        mom_aff = "1"
        dad_aff = "1"
        
        self.trio = self.create_family(sex, mom_aff, dad_aff)
        self.var = self.create_var(sex)
    
    def create_var(self, chrom='1', position='150', sex='F', child_geno='0/1'):
        ''' generate a test variant
        '''
        
        child = create_snv(sex, child_geno, chrom=chrom, pos=position)
        mom = create_snv("F", "0/0", chrom=chrom, pos=position)
        dad = create_snv("M", "0/0", chrom=chrom, pos=position)
        
        return TrioGenotypes(child.get_chrom(), child.get_position(),
            child, mom, dad)
    
    def create_family(self, child_gender, mom_aff, dad_aff):
        """ create a default family, with optional gender and parental statuses
        """
        
        fam = Family('test')
        fam.add_child('child', 'mother', 'father', child_gender, '2', 'child_vcf')
        fam.add_mother('mother', '0', '0', 'female', mom_aff, 'mother_vcf')
        fam.add_father('father', '0', '0', 'male', dad_aff, 'father_vcf')
        fam.set_child()
        
        return fam
    
    def test_getters(self):
        ''' test that the class getter methods work correctly
        '''
        
        # create a variant where the child, mother and father variables are filled in
        var = self.create_var(chrom='1', position='150', sex='F', child_geno='0/1')
        self.assertEqual(var.get_chrom(), '1')
        self.assertEqual(var.get_position(), 150)
        self.assertEqual(var.get_genes(), [['TEST']])
        self.assertEqual(var.get_range(), (150, 150))
        self.assertEqual(var.is_cnv(), False)
        self.assertEqual(var.get_inheritance_type(), 'autosomal')
        
        # construct a variant without values for the getter methods
        var = TrioGenotypes()
        self.assertEqual(var.get_chrom(), None)
        self.assertEqual(var.get_position(), None)
        self.assertEqual(var.get_genes(), None)
        self.assertEqual(var.get_range(), None)
        self.assertEqual(var.is_cnv(), None)
        self.assertEqual(var.get_inheritance_type(), None)
    
    def test_passes_de_novo_checks(self):
        """ test that passes_de_novo_checks() works correctly
        """
        
        # check that a default de novo variant passes
        self.assertTrue(self.var.passes_de_novo_checks(pp_filter=0.9))
        
        # check that vars fail without DENOVO-SNP or DENOVO-INDEL flags
        del self.var.child.info["DENOVO-SNP"]
        self.assertFalse(self.var.passes_de_novo_checks(pp_filter=0.9))
        
        # make sure that DENOVO-INDEL flag can pass the de novo filter
        self.var.child.info["DENOVO-INDEL"] = True
        self.assertTrue(self.var.passes_de_novo_checks(pp_filter=0.9))
        
        # check that de novos with low PP_DNM scores fail the de novo filter
        self.var.child.format["PP_DNM"] = 0.0099
        self.assertFalse(self.var.passes_de_novo_checks(pp_filter=0.9))
        
        # check that de novos with low PP_DNM scores pass the de novo filter, if
        # we are using a low PP_DNM threshold
        self.var.child.format["PP_DNM"] = 0.0099
        self.assertTrue(self.var.passes_de_novo_checks(pp_filter=0.0))
        
        # check that we don't fail a de novo if it lacks the PP_DNM annotation
        del self.var.child.format["PP_DNM"]
        self.assertTrue(self.var.passes_de_novo_checks(pp_filter=0.9))
    
    def test_passes_de_novo_checks_X_chrom(self):
        """ test that passes_de_novo_checks() works on the X chromosome
        """
        
        # check that a male X chrom de novo passes
        var = self.create_var(chrom='X', sex='M', child_geno='1/1')
        self.assertTrue(var.passes_de_novo_checks(pp_filter=0.9))
        
        # and change a field so that it would fail
        del var.child.info["DENOVO-SNP"]
        self.assertFalse(var.passes_de_novo_checks(pp_filter=0.9))
        
        # now check that a female X chrom de novo passes
        var = self.create_var(chrom='X', sex='F', child_geno='0/1')
        self.assertTrue(var.passes_de_novo_checks(pp_filter=0.9))
    
    def test_get_de_novo_genotype(self):
        """ check that get_de_novo_genotype() works correctly
        """
        
        var = self.create_var(chrom='1')
        self.assertEqual(var.get_inheritance_type(), 'autosomal')
        self.assertEqual(var.get_de_novo_genotype(), (1, 0, 0))
        
        var = self.create_var(chrom='X')
        self.assertEqual(var.get_inheritance_type(), 'XChrFemale')
        self.assertEqual(var.get_de_novo_genotype(), (1, 0, 0))
        
        # we double the alt count for males on the X, so a de novo genotype
        # differes from the other situations
        var = self.create_var(chrom='X', sex='M', child_geno='1/1')
        self.assertEqual(var.get_inheritance_type(), 'XChrMale')
        self.assertEqual(var.get_de_novo_genotype(), (2, 0, 0))
    
    def test_get_trio_genotype(self):
        """ test that get_trio_genotype() works correctly
        """
        
        # check that the defaul var gives the expected genotypes
        self.assertEqual(self.var.get_trio_genotype(), (1,0,0))
        
        # check that different genotypes still work
        self.var.mother.format["GT"] = "1/1"
        self.var.mother.set_genotype()
        self.assertEqual(self.var.get_trio_genotype(), (1,2,0))
        
        # check that probands only give NA genotypes for parents
        self.var.mother = None
        self.var.father = None
        self.assertEqual(self.var.get_trio_genotype(), (1, None, None))
    
    def test_chrom_to_int(self):
        """ test that chrom_to_int() works correctly
        """
        
        # check that an autosomal chrom works
        self.assertEqual(self.var.chrom_to_int("1"), 1)
        
        # check that an X chrom works
        self.assertEqual(self.var.chrom_to_int("X"), 23)
        
        # check that an X chrom works
        self.assertEqual(self.var.chrom_to_int("chrX"), 23)
        
        # check that a Y chrom works
        self.assertEqual(self.var.chrom_to_int("chrY"), 24)


if __name__ == '__main__':
    unittest.main()
