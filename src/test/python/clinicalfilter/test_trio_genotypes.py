""" unit testing of the TrioGenotypes class
"""

import unittest
import copy

from clinicalfilter.ped import Family
from clinicalfilter.ped import Person
from clinicalfilter.variant import Variant
from clinicalfilter.variant_cnv import CNV
from clinicalfilter.variant_snv import SNV
from clinicalfilter.vcf_info import VcfInfo
from clinicalfilter.trio_genotypes import TrioGenotypes


class TestTrioGenotypesPy(unittest.TestCase):
    """ test the Inheritance class
    """
    
    def setUp(self):
        """ define a family and variant, and start the Inheritance class
        """
        
        # generate a test family
        child_gender = "F"
        mom_aff = "1"
        dad_aff = "1"
        
        self.trio = self.create_family(child_gender, mom_aff, dad_aff)
        
        # generate a test variant
        child_var = self.create_snv(child_gender, "0/1")
        mom_var = self.create_snv("F", "0/0")
        dad_var = self.create_snv("M", "0/0")
        
        self.var = TrioGenotypes(child_var)
        self.var.add_mother_variant(mom_var)
        self.var.add_father_variant(dad_var)
    
    def create_snv(self, gender, genotype):
        """ create a default variant
        """
        
        chrom = "1"
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "G"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = SNV(chrom, pos, snp_id, ref, alt, filt)
        
        tags = {"gene": ["HGNC", "VGN", "GN"], "consequence": ["VCQ", "CQ"]}
        
        info = "HGNC=TEST;CQ=missense_variant;DENOVO-SNP;PP_DNM=0.99"
        keys = "GT:DP:TEAM29_FILTER:PP_DNM"
        values = genotype + ":50:PASS:0.99"
        
        var.add_info(info, tags)
        var.add_format(keys, values)
        var.set_gender(gender)
        var.set_genotype()
        
        return var
    
    def create_family(self, child_gender, mom_aff, dad_aff):
        """ create a default family, with optional gender and parental statuses
        """
        
        fam = Family("test")
        fam.add_child("child", "child_vcf", "2", child_gender)
        fam.add_mother("mother", "mother_vcf", mom_aff, "2")
        fam.add_father("father", "father_vcf", dad_aff, "1")
        fam.set_child()
        
        return fam
    
    def test_passes_de_novo_checks(self):
        """ test that passes_de_novo_checks() works correctly
        """
        
        # check that a default de novo variant passes
        self.assertTrue(self.var.passes_de_novo_checks())
        
        # check that TEAM29_FILTER != PASS fail
        self.var.child.format["TEAM29_FILTER"] = "fail"
        self.assertFalse(self.var.passes_de_novo_checks())
        
        # check that vars fail without the TEAM29_FILTER field
        del self.var.child.format["TEAM29_FILTER"]
        self.assertFalse(self.var.passes_de_novo_checks())
        
        # put the TEAM29_FILTER in, so later tests don't fail from its lack
        self.var.child.format["TEAM29_FILTER"] = "PASS"
        
        # check that vars fail without DENOVO-SNP or DENOVO-INDEL flags
        del self.var.child.info["DENOVO-SNP"]
        self.assertFalse(self.var.passes_de_novo_checks())
        
        # make sure that DENOVO-INDEL flag can pass the de novo filter
        self.var.child.info["DENOVO-INDEL"] = True
        self.assertTrue(self.var.passes_de_novo_checks())
        
        # check that de novos with low PP_DNM scores fail the de novo filter
        self.var.child.format["PP_DNM"] = 0.0099
        self.assertFalse(self.var.passes_de_novo_checks())
        
        # check that we don't fail a de novo if it lacks the PP_DNM annotation
        del self.var.child.format["PP_DNM"]
        self.assertTrue(self.var.passes_de_novo_checks())
    
    def test_passes_de_novo_checks_X_chrom(self):
        """ test that passes_de_novo_checks() works on the X chromosome
        """
        
        # check that a male X chrom de novo passes
        self.trio.child.gender = "M"
        self.var.inheritance_type = "XChrMale"
        self.var.child.format["GT"] = "1/1"
        self.var.child.set_genotype()
        self.assertTrue(self.var.passes_de_novo_checks())
        
        # and change a field so that it would fail
        del self.var.child.info["DENOVO-SNP"]
        self.assertFalse(self.var.passes_de_novo_checks())
        
        # and change the variant fom a male X de novo genotype
        self.var.child.format["GT"] = "1/0"
        self.var.child.set_genotype()
        self.assertTrue(self.var.passes_de_novo_checks())
        
        # now check that a female X chrom de novo passes
        self.trio.child.gender = "F"
        self.var.inheritance_type = "XChrFemale"
        self.var.child.set_genotype()
        self.var.child.info["DENOVO-SNP"] = True
        self.assertTrue(self.var.passes_de_novo_checks())
    
    def test_get_de_novo_genotype(self):
        """ check that get_de_novo_genotype() works correctly
        """
        
        self.var.inheritance_type = "autosomal"
        self.assertEqual(self.var.get_de_novo_genotype(), (1, 0, 0))
        
        self.var.inheritance_type = "XChrFemale"
        self.assertEqual(self.var.get_de_novo_genotype(), (1, 0, 0))
        
        # we double the alt count for males on the X, so a de novo genotype
        # differes from the other situations
        self.var.inheritance_type = "XChrMale"
        self.assertEqual(self.var.get_de_novo_genotype(), (2, 0, 0))
    
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
        del self.var.mother
        del self.var.father
        self.assertEqual(self.var.get_trio_genotype(), (1,"NA","NA"))
    
    def test_convert_chrom_to_int(self):
        """ test that convert_chrom_to_int() works correctly
        """
        
        # check that an autosomal chrom works
        self.assertEqual(self.var.convert_chrom_to_int("1"), 1) 
        
        # check that an X chrom works
        self.assertEqual(self.var.convert_chrom_to_int("X"), 23)
        
        # check that an X chrom works
        self.assertEqual(self.var.convert_chrom_to_int("chrX"), 23)
        
        # check that a Y chrom works
        self.assertEqual(self.var.convert_chrom_to_int("chrY"), 24)


if __name__ == '__main__':
    unittest.main()


