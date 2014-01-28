""" unit testing of the SNV class
"""

import unittest
from clinicalfilter.variant_snv import SNV

class TestVariantSnvPy(unittest.TestCase):
    """
    """
    
    def setUp(self):
        """ define a default VcfInfo object
        """
        
        chrom = "1"
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "G"
        qual = "50"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        self.var = SNV(chrom, pos, snp_id, ref, alt, qual, filt)
        
        tags = {"gene": ["HGNC", "VGN", "GN"], "consequence": ["VCQ", "CQ"]}
        
        info = "HGNC=ATRX;CQ=missense_variant;random_tag"
        self.format_keys = "GT:DP"
        self.sample_values = "0/1:50"
        
        self.var.add_info(info, tags)
    
    def test_get_key(self):
        """ tests that get_key() operates correctly
        """
        
        # make sure the chrom and position are correct
        self.var.chrom = "1"
        self.var.position = "15000000"
        
        self.assertEqual(self.var.get_key(), ("1", "15000000"))
        
        # and make sure the chrom and position are correct if we change them
        self.var.chrom = "22"
        self.var.position = "123456789"
        self.assertEqual(self.var.get_key(), ("22", "123456789"))
    
    def test_convert_genotype(self):
        """ test that genotypes convert from two char to single char
        """
        
        genotypes = [("0/0", 0), ("0/1", 1), ("1/0", 1), ("1/1", 2), \
            ("1/2", 1), ("2/1", 1), ("0/2", 1), ("2/0", 1), ("2/2", 2)]
        
        # run thorugh all the legit genotype codes
        for geno in genotypes:
            genotype = geno[0]
            result = geno[1]
            self.assertEqual(self.var.convert_genotype(genotype), result)
         
        # Raise error when converting single character genotype
        with self.assertRaises(ValueError):
            self.var.convert_genotype("0")
          
        # raise error when converting unknown genotype
        with self.assertRaises(KeyError):
            self.var.convert_genotype("a/a")
            
        # also include other genotype format posibilities. None of these are
        # used, but since they aren't explicitly forbidden, make sure they work
        
        # check two character strings
        self.assertEqual(self.var.convert_genotype("00"), 0)
        self.assertEqual(self.var.convert_genotype("01"), 1)
        
        # check > three character strings
        self.assertEqual(self.var.convert_genotype("00001"), 1)
        self.assertEqual(self.var.convert_genotype("1000001"), 2)
    
    def test_set_default_genotype(self):
        """ test that set_default_genotype() operates correctly on the autosomes
        """
        
        self.var.set_gender("male")
        self.var.set_default_genotype()
        self.assertEqual(self.var.get_genotype(), 0)
    
    def test_set_genotype_autosomal(self):
        """ test that set_genotype() operates correctly
        """
        
        self.var.add_format(self.format_keys, self.sample_values)
        self.var.set_gender("male")
        
        genotypes = [("0/0", 0), ("0/1", 1), ("1/1", 2)]
        
        for geno in genotypes:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.get_genotype(), result)
        
        # remove the format attribute, so we can raise an error
        del self.var.format
        with self.assertRaises(ValueError):
            self.var.set_genotype()
    
    def test_set_genotype_allosomal_male(self):
        """ test that set_genotype() operates correctly for the male X chrom
        """
        
        self.var.add_format(self.format_keys, self.sample_values)
        self.var.chrom = "X"
        self.var.set_gender("male")
        
        genotypes = [("0/0", 0), ("1/1", 2)]
        
        for geno in genotypes:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.get_genotype(), result)
        
        # check that we raise an error for X chrom hets
        genotypes = ["0/1", "1/0"]
        for genotype in genotypes:
            self.var.format["GT"] = genotype
            with self.assertRaises(ValueError):
                self.var.set_genotype()
    
    def test_set_genotype_allosomal_female(self):
        """ test that set_genotype() operates correctly for the female X chrom
        """
        
        self.var.add_format(self.format_keys, self.sample_values)
        self.var.chrom = "X"
        self.var.set_gender("female")
        
        genotypes = [("0/0", 0), ("0/1", 1), ("1/1", 2)]
        
        for geno in genotypes:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.get_genotype(), result)
    
    def test_is_het_autosomal(self):
        """ tests that is_het() operates correctly for automsal chromosomes
        """
        
        self.var.add_format(self.format_keys, self.sample_values)
        self.var.set_gender("male")
        
        het = [("0/0", False), ("0/1", True), ("1/1", False)]
        
        for geno in het:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.is_het(), result)
     
    def test_is_hom_alt_autosomal(self):
        """ tests that is_hom_alt() operates correctly for automsal chromosomes
        """
        
        self.var.add_format(self.format_keys, self.sample_values)
        self.var.set_gender("male")
        
        hom_alt = [("0/0", False), ("0/1", False), ("1/1", True)]
        
        for geno in hom_alt:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.is_hom_alt(), result)
            
    def test_is_hom_ref_autosomal(self):
        """ tests that is_hom_ref() operates correctly for automsal chromosomes
        """
        
        self.var.add_format(self.format_keys, self.sample_values)
        self.var.set_gender("male")
        
        hom_ref = [("0/0", True), ("0/1", False), ("1/1", False)]
        
        for geno in hom_ref:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.is_hom_ref(), result)
        
    def test_is_not_ref_autosomal(self):
        """ tests that is_not_ref() operates correctly for automsal chromosomes
        """
        
        self.var.add_format(self.format_keys, self.sample_values)
        self.var.set_gender("male")
        
        not_ref = [("0/0", False), ("0/1", True), ("1/1", True)]
        
        for geno in not_ref:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.is_not_ref(), result)
    
    def test_is_not_alt_autosomal(self):
        """ tests that is_not_ref() operates correctly for automsal chromosomes
        """
        
        self.var.add_format(self.format_keys, self.sample_values)
        self.var.set_gender("male")
        
        not_alt = [("0/0", True), ("0/1", True), ("1/1", False)]
        
        for geno in not_alt:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.is_not_alt(), result)
    
    

unittest.main()

