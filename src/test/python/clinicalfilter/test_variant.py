""" unit testing of the Variant class
"""

import unittest
from clinicalfilter.variant import Variant


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
        qual = "50"
        filt = "PASS"
        
        self.var = Variant(chrom, pos, snp_id, ref, alt, qual, filt)
    
    def test_set_gender_unknown(self):
        """ tests set_gender(), and its implications on the inheritance type
        """
        
        # raise error for unknown gender
        with self.assertRaises(ValueError):
            self.var.set_gender("unknown")
    
    def test_set_gender_autosomal(self):
        """test gender on autosomal chroms
        """
        
        # check all the legitimate gender codes
        gender_codes = ["1", "m", "M", "male", "2", "f", "F", "female"]
        for gender in gender_codes:
            self.var.set_gender(gender)
            self.assertEqual(self.var.get_inheritance_type(), "autosomal")
    
    def test_set_gender_allosomal(self):
        """ test gender on the X chromosome
        """
        
        # check the allosomal chroms
        self.var.chrom = "X"
        self.var.set_gender("male")
        self.assertEqual(self.var.get_inheritance_type(), "XChrMale")
        
        self.var.set_gender("female")
        self.assertEqual(self.var.get_inheritance_type(), "XChrFemale")
    
    def test_set_gender_pseudoautosomal(self):
        """ test gender on a pseudoautosomal region
        """
        # now check a variant in the pseudoautosomal regions
        var = self.var
        self.var.chrom = "X"
        self.var.position = "2699510"
        self.var.set_gender("male")
        self.assertEqual(self.var.get_inheritance_type(), "autosomal")
        
        self.var.set_gender("female")
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
        self.assertEqual(self.var.get_position(), "15000000")
        
        # check that we can change the chrom and it still returns correctly
        self.var.position = "123456789"
        self.assertEqual(self.var.get_position(), "123456789")
    
    def test_get_mutation_id(self):
        """ test that the mutation ID is parsed correctly
        """
        var = self.var
        
        # check a null ID value
        var.id = "."
        var.set_mutation_id()
        self.assertEqual(var.get_mutation_id(), "NA")
        
        # check a rsID value
        var.id = "rs12546"
        var.set_mutation_id()
        self.assertEqual(var.get_mutation_id(), "NA")
        
        # check a rsID and a mutation ID
        var.id = "rs12546&CM0001"
        var.set_mutation_id()
        self.assertEqual(var.get_mutation_id(), "CM0001")
        
        # check multiple mutation IDs
        var.id = "CM0001&CM0002"
        var.set_mutation_id()
        self.assertEqual(var.get_mutation_id(), "CM0001,CM0002")
    
    def test_get_vcf_line(self):
        """ tests that the vcf line string returns correctly
        """
        
        vcf_line = ["1", "15000000", ".", "A", "G", "50", "PASS", "AB=0.41;AC=1;AN=2", "GT:gatk_PL:GQ", "0/1:736,0,356:99"]
        self.var.add_vcf_line(vcf_line)
        self.assertEqual(self.var.get_vcf_line(), vcf_line)
    
    # TODO: check add_format


# unittest.main()

