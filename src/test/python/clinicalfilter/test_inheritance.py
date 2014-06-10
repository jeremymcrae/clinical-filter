""" unit testing of the Inheritance class
"""

import unittest
import sys

from clinicalfilter.ped import Family
from clinicalfilter.ped import Person
from clinicalfilter.variant import Variant
from clinicalfilter.variant_snv import SNV
from clinicalfilter.variant_cnv import CNV
from clinicalfilter.inheritance import Autosomal
from clinicalfilter.inheritance import Allosomal
from clinicalfilter.vcf_info import VcfInfo
from clinicalfilter.trio_genotypes import TrioGenotypes

IS_PYTHON2 = sys.version_info[0] == 2
IS_PYTHON3 = sys.version_info[0] == 3


class TestInheritancePy(unittest.TestCase):
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
        
        # generate list of variants
        self.variants = [self.create_variant(child_gender)]
        self.variants.append(self.create_variant(child_gender))
        
        # make sure we've got known genes data
        self.known_genes = {"TEST": {"inheritance": ["Monoallelic"], "confirmed_status": ["Confirmed DD Gene"]}}
        gene_inh = self.known_genes[self.variants[0].get_gene()]["inheritance"]
        
        self.inh = Autosomal(self.variants, self.trio, gene_inh)
    
    def create_snv(self, gender, genotype, chrom, pos):
        """ create a default variant
        """
        
        snp_id = "."
        ref = "A"
        alt = "G"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = SNV(chrom, pos, snp_id, ref, alt, filt)
        
        tags = {"gene": ["HGNC", "VGN", "GN"], "consequence": ["VCQ", "CQ"]}
        
        info = "HGNC=TEST;CQ=missense_variant;random_tag"
        format_keys = "GT:DP"
        sample_values = genotype + ":50"
        
        var.add_info(info, tags)
        var.add_format(format_keys, sample_values)
        var.set_gender(gender)
        var.set_genotype()
        
        return var
    
    def create_cnv(self, gender, inh, chrom, pos):
        """ create a default variant
        """
        
        snp_id = "."
        ref = "A"
        alt = "<DUP>"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = CNV(chrom, pos, snp_id, ref, alt, filt)
        
        tags = {"gene": ["HGNC", "VGN", "GN"], "consequence": ["VCQ", "CQ"]}
        
        info = "HGNC=TEST;HGNC_ALL=TEST;END=16000000;SVLEN=5000"
        format_keys = "INHERITANCE:DP"
        sample_values = inh + ":50"
        
        var.add_info(info, tags)
        var.add_format(format_keys, sample_values)
        var.set_gender(gender)
        var.set_genotype()
        
        return var
    
    def create_variant(self, child_gender, chrom="1", position="15000000"):
        """ creates a TrioGenotypes variant
        """
        
        # generate a test variant
        child_var = self.create_snv(child_gender, "0/1", chrom, position)
        mom_var = self.create_snv("F", "0/0", chrom, position)
        dad_var = self.create_snv("M", "0/0", chrom, position)
        
        var = TrioGenotypes(child_var)
        var.add_mother_variant(mom_var)
        var.add_father_variant(dad_var)
        
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
    
    def test_check_inheritance_mode_matches_gene_mode(self):
        """ test that check_inheritance_mode_matches_gene_mode() works correctly
        """
        
        # check that the default inheritance types have been set up correctly
        self.assertEqual(self.inh.inheritance_modes, {"Monoallelic", "Biallelic", "Both"})
        
        # make sure that the default var and gene inheritance work
        self.assertTrue(self.inh.check_inheritance_mode_matches_gene_mode())
        
        # check that no gene inheritance overlap fails
        self.inh.gene_inheritance = {"Mosaic"}
        self.inh.inheritance_modes = {"Monoallelic", "Biallelic", "Both"}
        self.assertFalse(self.inh.check_inheritance_mode_matches_gene_mode())
        
        # check that a single inheritance type still works
        self.inh.gene_inheritance = {"Monoallelic"}
        self.assertTrue(self.inh.check_inheritance_mode_matches_gene_mode())
        
        # check that multiple inheritance types for a gene still work
        self.inh.gene_inheritance = {"Monoallelic", "Biallelic"}
        self.assertTrue(self.inh.check_inheritance_mode_matches_gene_mode())
        
        # check that extra inheritance modes are included still work
        self.inh.gene_inheritance = {"Monoallelic", "Biallelic", "Mosaic"}
        self.assertTrue(self.inh.check_inheritance_mode_matches_gene_mode())
    
    def test_set_trio_genotypes(self):
        """ test that set_trio_genotypes() works correctly
        """
        
        # set the genotypes using the default variant
        var = self.variants[0]
        self.inh.set_trio_genotypes(var)
        
        # the genotypes for the inh object should match the vars genotypes
        self.assertEqual(self.inh.child, var.child)
        self.assertEqual(self.inh.mom, var.mother)
        self.assertEqual(self.inh.dad, var.father)
        
        # now remove the parents before re-setting the genotypes
        del var.mother
        del var.father
        self.inh.trio.father = None
        self.inh.trio.mother = None
        self.inh.set_trio_genotypes(var)
        
        # the child should match the vars genotypes, but the parent's 
        # genotypes should be None
        self.assertEqual(self.inh.child, var.child)
        self.assertIsNone(self.inh.mom)
        self.assertIsNone(self.inh.dad)
    
    def test_add_variant_to_appropriate_list(self):
        """ test that add_variant_to_appropriate_list() works correctly
        """
        
        var = self.variants[0]
        inheritance = "Monoallelic"
        check = "compound_het"
        
        # check that compound_het vars are only added to the compound_het list
        self.inh.compound_hets = []
        self.inh.candidates = []
        self.inh.add_variant_to_appropriate_list(var, check, inheritance)
        self.assertEqual(self.inh.candidates, [])
        self.assertEqual(self.inh.compound_hets, [(var, check, inheritance)])
        
        # check that single_variant vars are only added to the candidates list
        self.inh.compound_hets = []
        self.inh.candidates = []
        check = "single_variant"
        self.inh.add_variant_to_appropriate_list(var, check, inheritance)
        self.assertEqual(self.inh.candidates, [(var, check, inheritance)])
        self.assertEqual(self.inh.compound_hets, [])
        
        # check that other vars aren't added either list
        self.inh.compound_hets = []
        self.inh.candidates = []
        check = "nothing"
        self.inh.add_variant_to_appropriate_list(var, check, inheritance)
        self.assertEqual(self.inh.candidates, [])
        self.assertEqual(self.inh.compound_hets, [])
    
    def test_check_if_any_variant_is_cnv(self):
        """ test if check_if_any_variant_is_cnv() works correctly
        """
        
        # generate a test variant
        chrom = "1"
        position = "60000"
        child_var = self.create_cnv("F", "unknown", chrom, position)
        mom_var = self.create_cnv("F", "unknown", chrom, position)
        dad_var = self.create_cnv("M", "unknown", chrom, position)
        
        cnv_var = TrioGenotypes(child_var)
        cnv_var.add_mother_variant(mom_var)
        cnv_var.add_father_variant(dad_var)
        
        # check that all variants=SNV returns False
        self.assertFalse(self.inh.check_if_any_variant_is_cnv())
        
        # add a CNV to the variants, then check that we find a CNV
        self.inh.variants.append(cnv_var)
        self.assertTrue(self.inh.check_if_any_variant_is_cnv())
    
    def set_compound_het_var(self, var, geno, compound_type):
        """ convenience function to set the trio genotypes for a variant
        """
        
        genos = {"0": "0/0", "1": "0/1", "2": "1/1"}
        
        # convert the geno codes to allele codes
        child = genos[geno[0]]
        mom = genos[geno[1]]
        dad = genos[geno[2]]
        
        # set the genotype field for each individual
        var.child.format["GT"] = child
        var.mother.format["GT"] = mom
        var.father.format["GT"] = dad
        
        # and set th genotype for each individual
        var.child.set_genotype()
        var.mother.set_genotype()
        var.father.set_genotype()
        
        # set the trio genotypes for the inheritance object
        return (var, compound_type, "Biallelic")
    
    def test_check_compound_hets_autosomal(self):
        """ test that check_compound_hets() works correctly for autosomal vars
        """
        
        # set some variants, so we can alter them later
        var1 = self.create_variant("F", chrom="1", position="15000000")
        var2 = self.create_variant("F", chrom="1", position="16000000")
        var3 = self.create_variant("F", chrom="1", position="17000000")
        
        # set the inheritance type, the compound het type ("compound_het" 
        # for autosomal variants, and start autosomal inheritance)
        inh = "Biallelic"
        compound = "compound_het"
        self.inh = Autosomal([var1, var2, var3], self.trio, inh)
        
        variants = ["", ""]
        
        # check the expected "110, 101" combo passes
        variants[0] = self.set_compound_het_var(var1, "110", compound)
        variants[1] = self.set_compound_het_var(var2, "101", compound)
        if IS_PYTHON3:
            self.assertCountEqual(self.inh.check_compound_hets(variants), variants)
        elif IS_PYTHON2:
            self.assertEqual(sorted(self.inh.check_compound_hets(variants)), sorted(variants))
        
        # check that "110, 110" combo fails
        variants[0] = self.set_compound_het_var(var1, "110", compound)
        variants[1] = self.set_compound_het_var(var2, "110", compound)
        self.assertEqual(self.inh.check_compound_hets(variants), [])
        
        # check that "101, 101" combo fails
        variants[0] = self.set_compound_het_var(var1, "101", compound)
        variants[1] = self.set_compound_het_var(var2, "101", compound)
        self.assertEqual(self.inh.check_compound_hets(variants), [])
        
        # check that > 2 valid compound hets passes all variants
        variants = ["", "", ""]
        variants[0] = self.set_compound_het_var(var1, "110", compound)
        variants[1] = self.set_compound_het_var(var2, "101", compound)
        variants[2] = self.set_compound_het_var(var3, "110", compound)
        if IS_PYTHON3:
            self.assertCountEqual(self.inh.check_compound_hets(variants), variants)
        elif IS_PYTHON2:
            self.assertEqual(sorted(self.inh.check_compound_hets(variants)), sorted(variants))
        
        # check that a single var fails to give compound hets
        single_var = variants[:1]
        self.assertEqual(self.inh.check_compound_hets(single_var), [])
        
        # check that zero length list gives no compound hets
        no_vars = []
        self.assertEqual(self.inh.check_compound_hets(no_vars), [])
        
        # check that de novo containing "110, 100" combos give compound hets
        variants = ["", ""]
        variants[0] = self.set_compound_het_var(var1, "110", compound)
        variants[1] = self.set_compound_het_var(var2, "100", compound)
        if IS_PYTHON3:
            self.assertCountEqual(self.inh.check_compound_hets(variants), variants)
        elif IS_PYTHON2:
            self.assertEqual(sorted(self.inh.check_compound_hets(variants)), sorted(variants))
        
        # check that de novo "100, 100" combos give compound hets
        variants[0] = self.set_compound_het_var(var1, "100", compound)
        variants[1] = self.set_compound_het_var(var2, "100", compound)
        if IS_PYTHON3:
            self.assertCountEqual(self.inh.check_compound_hets(variants), variants)
        elif IS_PYTHON2:
            self.assertEqual(sorted(self.inh.check_compound_hets(variants)), sorted(variants))
        
        # check that "111, 111" combos require affected parents
        variants[0] = self.set_compound_het_var(var1, "111", compound)
        variants[1] = self.set_compound_het_var(var2, "111", compound)
        self.assertEqual(self.inh.check_compound_hets(variants), [])
        
        # check "111, 111" combo with a single affected parent
        self.inh.mother_affected = True
        self.assertEqual(self.inh.check_compound_hets(variants), [])
        
        # check "111, 111" combo with both parents affected
        self.inh.father_affected = True
        if IS_PYTHON3:
            self.assertCountEqual(self.inh.check_compound_hets(variants), variants)
        elif IS_PYTHON2:
            self.assertEqual(sorted(self.inh.check_compound_hets(variants)), sorted(variants))
        
        # check that without parents, all variants are included, even if they 
        # wouldn't pass normally
        self.inh.trio.mother = None
        self.inh.trio.father = None
        variants[0] = self.set_compound_het_var(var1, "101", compound)
        variants[1] = self.set_compound_het_var(var2, "101", compound)
        if IS_PYTHON3:
            self.assertCountEqual(self.inh.check_compound_hets(variants), variants)
        elif IS_PYTHON2:
            self.assertEqual(sorted(self.inh.check_compound_hets(variants)), sorted(variants))
        
    def test_check_compound_hets_allosomal(self):
        """ test that check_compound_hets() works correctly for allosomal vars
        """
        
        # set some X chrom variants, so we can alter them later
        var1 = self.create_variant("F", chrom="X", position="15000000")
        var2 = self.create_variant("F", chrom="X", position="16000000")
        
        # set the inheritance type, the compound het type ("hemizygous" 
        # for allosomal variants, and start allosomal inheritance)
        inh = "Hemizygous"
        compound = "hemizygous"
        self.inh = Allosomal([var1, var2], self.trio, inh)
        
        variants = ["", ""]
        
        # check that de novo containing "110, 100" combos pass
        variants[0] = self.set_compound_het_var(var1, "110", compound)
        variants[1] = self.set_compound_het_var(var2, "100", compound)
        if IS_PYTHON3:
            self.assertCountEqual(self.inh.check_compound_hets(variants), variants)
        elif IS_PYTHON2:
            self.assertEqual(sorted(self.inh.check_compound_hets(variants)), sorted(variants))
        
        # check that "110, 102" combo fails if the father is unaffected
        variants[0] = self.set_compound_het_var(var1, "110", compound)
        variants[1] = self.set_compound_het_var(var2, "102", compound)
        self.assertEqual(self.inh.check_compound_hets(variants), [])
        
        # check that "110, 102" combo passes if the father is affected
        self.inh.father_affected = True
        if IS_PYTHON3:
            self.assertCountEqual(self.inh.check_compound_hets(variants), variants)
        elif IS_PYTHON2:
            self.assertEqual(sorted(self.inh.check_compound_hets(variants)), sorted(variants))
        
        # make sure we can't set the father as het on the X chrom
        with self.assertRaises(ValueError):
            self.set_compound_het_var(var2, "101", compound)


if __name__ == '__main__':
    unittest.main()

