""" unit testing of the Inheritance class
"""

import unittest

from clinicalfilter.ped import Family
from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.inheritance import Autosomal
from clinicalfilter.inheritance import Allosomal
from clinicalfilter.trio_genotypes import TrioGenotypes


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
        self.known_genes = {"TEST": {"inh": ["Monoallelic"], "confirmed_status": ["Confirmed DD Gene"]}}
        
        self.inh = Autosomal(self.variants, self.trio, self.known_genes, "TEST")
    
    def create_snv(self, gender, genotype, chrom, pos, cq=None):
        """ create a default variant
        """
        
        snp_id = "."
        ref = "A"
        alt = "G"
        filt = "PASS"
        
        if cq is None:
            cq = "missense_variant"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = SNV(chrom, pos, snp_id, ref, alt, filt)
        
        info = "HGNC=TEST;CQ={};random_tag".format(cq)
        format_keys = "GT:DP"
        sample_values = genotype + ":50"
        
        var.add_info(info)
        var.add_format(format_keys, sample_values)
        var.set_gender(gender)
        var.set_genotype()
        
        return var
    
    def create_cnv(self, gender, inh, chrom, pos, cq=None):
        """ create a default variant
        """
        
        snp_id = "."
        ref = "A"
        alt = "<DEL>"
        filt = "PASS"
        
        if cq is None:
            cq = "transcript_ablation"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = CNV(chrom, pos, snp_id, ref, alt, filt)
        
        info = "CQ={};HGNC=TEST;HGNC_ALL=TEST;END=16000000;SVLEN=5000".format(cq)
        format_keys = "INHERITANCE:DP:CIFER_INHERITANCE"
        sample_values = "{0}:50:{0}".format(inh)
        
        var.add_info(info)
        var.add_format(format_keys, sample_values)
        var.set_gender(gender)
        var.set_genotype()
        
        return var
    
    def create_variant(self, child_gender, chrom="1", position="15000000", cq=None):
        """ creates a TrioGenotypes variant
        """
        
        # generate a test variant
        try:
            child_var = self.create_snv(child_gender, "0/1", chrom, position, cq)
        except ValueError:
            child_var = self.create_snv(child_gender, "1/1", chrom, position, cq)
        
        mom_var = self.create_snv("F", "0/0", chrom, position, cq)
        dad_var = self.create_snv("M", "0/0", chrom, position, cq)
        
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
        self.assertEqual(self.inh.compound_hets, [(var, (check,), (inheritance,))])
        
        # check that single_variant vars are only added to the candidates list
        self.inh.compound_hets = []
        self.inh.candidates = []
        check = "single_variant"
        self.inh.add_variant_to_appropriate_list(var, check, inheritance)
        self.assertEqual(self.inh.candidates, [(var, [check], [inheritance])])
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
    
    def set_compound_het_var(self, var, geno):
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
        return var
    
    def test_check_compound_hets(self):
        """ test that check_compound_hets() works correctly for autosomal vars
        """
        
        # set some variants, so we can alter them later
        var1 = self.create_variant("F", chrom="1", position="15000000", cq="stop_gained")
        var2 = self.create_variant("F", chrom="1", position="16000000", cq="stop_gained")
        var3 = self.create_variant("F", chrom="1", position="17000000", cq="stop_gained")
        
        # set the inheritance type, the compound het type ("compound_het"
        # for autosomal variants, and start autosomal inheritance)
        # known_genes = "Biallelic"
        known_genes = {"TEST": {"inh": ["Biallelic"], "confirmed_status": ["Confirmed DD Gene"]}}
        self.inh = Autosomal([var1, var2, var3], self.trio, known_genes, "TEST")
        
        variants = [(), ()]
        
        # check the expected "110, 101" combo passes
        variants[0] = (self.set_compound_het_var(var1, "110"),)
        variants[1] = (self.set_compound_het_var(var2, "101"),)
        self.assertEqual(sorted(self.inh.check_compound_hets(variants)), sorted(variants))
        
        # check that > 2 valid compound hets passes all variants
        variants = [(), (), ()]
        variants[0] = (self.set_compound_het_var(var1, "110"),)
        variants[1] = (self.set_compound_het_var(var2, "101"),)
        variants[2] = (self.set_compound_het_var(var3, "110"),)
        self.assertEqual(sorted(self.inh.check_compound_hets(variants)), sorted(variants))
        
        # check that a single var fails to give compound hets
        single_var = variants[:1]
        self.assertEqual(self.inh.check_compound_hets(single_var), [])
        
        # check that zero length list gives no compound hets
        no_vars = []
        self.assertEqual(self.inh.check_compound_hets(no_vars), [])
    
    def test_is_compound_pair_identical_variants(self):
        """ check that is_compound_pair() excludes compound pairs where the
        members are identical
        """
        
        # set some variants, so we can alter them later
        var1 = self.create_variant("F", chrom="1", position="150", cq="stop_gained")
        var2 = self.create_variant("F", chrom="1", position="160", cq="stop_gained")
        
        var1 = self.set_compound_het_var(var1, "110")
        var2 = self.set_compound_het_var(var2, "101")
        
        # don't include pairs where the first and the second variant are identical
        self.assertFalse(self.inh.is_compound_pair(var1, var1))
        
        # make sure it works normally
        self.assertTrue(self.inh.is_compound_pair(var1, var2))
    
    def test_is_compound_pair_both_missense_with_parents(self):
        """check that is_compound_pair() excludes pairs where both are missense
        """
        
        # set some variants, so we can alter them later
        var1 = self.create_variant("F", chrom="1", position="150", cq="missense_variant")
        var2 = self.create_variant("F", chrom="1", position="160", cq="missense_variant")
        var3 = self.create_variant("F", chrom="1", position="160", cq="inframe_deletion")
        var4 = self.create_variant("F", chrom="1", position="160", cq="stop_gained")
        
        var1 = self.set_compound_het_var(var1, "110")
        var2 = self.set_compound_het_var(var2, "101")
        var3 = self.set_compound_het_var(var3, "101")
        var4 = self.set_compound_het_var(var4, "101")
        
        # dont exclude pairs where both members are not loss-of-function if the
        # proband has parents
        self.assertTrue(self.inh.is_compound_pair(var1, var2))
        self.assertTrue(self.inh.is_compound_pair(var1, var3))
        
        # make sure it works normally
        self.assertTrue(self.inh.is_compound_pair(var1, var4))
    
    def test_is_compound_pair_both_missense_without_parents(self):
        """check that is_compound_pair() excludes pairs where both are missense
        """
        
        # set some variants, so we can alter them later
        var1 = self.create_variant("F", chrom="1", position="150", cq="missense_variant")
        var2 = self.create_variant("F", chrom="1", position="160", cq="missense_variant")
        var3 = self.create_variant("F", chrom="1", position="160", cq="inframe_deletion")
        var4 = self.create_variant("F", chrom="1", position="160", cq="stop_gained")
        
        var1 = self.set_compound_het_var(var1, "110")
        var2 = self.set_compound_het_var(var2, "101")
        var3 = self.set_compound_het_var(var3, "101")
        var4 = self.set_compound_het_var(var4, "101")
        
        # drop the parents
        self.inh.trio.father = None
        self.inh.trio.mother = None
        
        # exclude pairs where both members are not loss-of-function
        self.assertFalse(self.inh.is_compound_pair(var1, var2))
        self.assertFalse(self.inh.is_compound_pair(var1, var3))
        
        # make sure it works if one variant is loss-of-function
        self.assertTrue(self.inh.is_compound_pair(var1, var4))
    
    def test_is_compound_pair_unknown_gene(self):
        """check that is_compound_pair() excludes pairs for unknown genes
        """
        
        # set some variants, so we can alter them later
        var1 = self.create_variant("F", chrom="1", position="150", cq="stop_gained")
        var2 = self.create_variant("F", chrom="1", position="160", cq="stop_gained")
        
        var1 = self.set_compound_het_var(var1, "110")
        var2 = self.set_compound_het_var(var2, "101")
        
        var1.child.genes = ["."]
        var2.child.genes = ["."]
        
        # exclude pairs where both members are not loss-of-function
        self.assertFalse(self.inh.is_compound_pair(var1, var2))
    
    def test_is_compound_pair_cnv_paternal(self):
        """ check that is_compound_pair() includes pairs with CNVs
        """
        
        # generate a test variant
        chrom = "1"
        position = "60000"
        child_var = self.create_cnv("F", "paternal", chrom, position)
        mom_var = self.create_cnv("F", "unknown", chrom, position)
        dad_var = self.create_cnv("M", "unknown", chrom, position)
        
        cnv = TrioGenotypes(child_var)
        cnv.add_mother_variant(mom_var)
        cnv.add_father_variant(dad_var)
        
        # set some variants, so we can alter them later
        snv = self.create_variant("F", chrom="1", position="150", cq="stop_gained")
        snv = self.set_compound_het_var(snv, "110")
        
        # check that these variants are compound hets, no matter which order
        # they are given as.
        self.assertTrue(self.inh.is_compound_pair(cnv, snv))
        self.assertTrue(self.inh.is_compound_pair(snv, cnv))
        
        # check that if the SNV is inherited from the same parent as the CNV,
        # then the pair isn't a compound het.
        snv = self.set_compound_het_var(snv, "101")
        self.assertFalse(self.inh.is_compound_pair(cnv, snv))
    
    def test_is_compound_pair_cnv_maternal(self):
        """ check that is_compound_pair() includes pairs with CNVs
        """
        
        # generate a test variant
        chrom = "1"
        position = "60000"
        child_var = self.create_cnv("F", "maternal", chrom, position)
        mom_var = self.create_cnv("F", "unknown", chrom, position)
        dad_var = self.create_cnv("M", "unknown", chrom, position)
        
        cnv = TrioGenotypes(child_var)
        cnv.add_mother_variant(mom_var)
        cnv.add_father_variant(dad_var)
        
        # set some variants, so we can alter them later
        snv = self.create_variant("F", chrom="1", position="150", cq="stop_gained")
        snv = self.set_compound_het_var(snv, "101")
        
        # check that these variants are compound hets, no matter which order
        # they are given as.
        self.assertTrue(self.inh.is_compound_pair(cnv, snv))
        
        # check that if the SNV is inherited from the same parent as the CNV,
        # then the pair isn't a compound het.
        snv = self.set_compound_het_var(snv, "110")
        self.assertFalse(self.inh.is_compound_pair(cnv, snv))
    
    def test_is_compound_pair_proband_only(self):
        """ check that is_compound_pair() includes proband-only pairs
        """
        
        fam = Family("test")
        fam.add_child("child", "child_vcf", "2", "F")
        fam.set_child()
        
        # set some variants, so we can alter them later
        var1 = self.create_variant("F", chrom="1", position="150", cq="stop_gained")
        var2 = self.create_variant("F", chrom="1", position="160", cq="stop_gained")
        
        inh = Autosomal([var1, var2], fam, self.known_genes, "TEST")
        
        # check that a proband-only passes, regardless of the parental genotypes
        self.assertTrue(inh.is_compound_pair(var1, var2))
    
    def test_is_compound_pair_allosomal(self):
        """ check that is_compound_pair() works when the father is affected
        """
        
        # set some variants, so we can alter them later
        var1 = self.create_variant("M", chrom="X", position="150", cq="stop_gained")
        var2 = self.create_variant("M", chrom="X", position="160", cq="stop_gained")
        
        var1 = self.set_compound_het_var(var1, "210")
        var2 = self.set_compound_het_var(var2, "202")
        
        # check when the father is unaffected
        self.assertFalse(self.inh.is_compound_pair(var1, var2))
        
        # now check when the father is affected
        self.inh.father_affected = True
        self.assertTrue(self.inh.is_compound_pair(var1, var2))
        
        # make sure we can't set the father as het on the X chrom
        with self.assertRaises(ValueError):
            self.set_compound_het_var(var2, "201")
    
    def test_is_compound_pair_genotype_combinations(self):
        """ check the various genotype combinations for a compound het
        """
        
        # set some variants, so we can alter them later
        var1 = self.create_variant("M", chrom="1", position="150", cq="stop_gained")
        var2 = self.create_variant("M", chrom="1", position="160", cq="stop_gained")
        
        var1 = self.set_compound_het_var(var1, "110")
        var2 = self.set_compound_het_var(var2, "101")
        
        # check that the expected scenario passes
        self.assertTrue(self.inh.is_compound_pair(var1, var2))
        
        # check that compound hets with de novos fail
        var1 = self.set_compound_het_var(var1, "100")
        var2 = self.set_compound_het_var(var2, "101")
        self.assertFalse(self.inh.is_compound_pair(var1, var2))
        
        # check that compound hets have to be transmitted from both parents
        var1 = self.set_compound_het_var(var1, "101")
        var2 = self.set_compound_het_var(var2, "101")
        self.assertFalse(self.inh.is_compound_pair(var1, var2))
        
        # check that compound hets have to be transmitted from both parents
        var1 = self.set_compound_het_var(var1, "110")
        var2 = self.set_compound_het_var(var2, "110")
        self.assertFalse(self.inh.is_compound_pair(var1, var2))
        
        # make sure that only compound hets in trans pass. We exclude
        var1 = self.set_compound_het_var(var1, "111")
        var2 = self.set_compound_het_var(var2, "101")
        self.assertFalse(self.inh.is_compound_pair(var1, var2))
        
        var1 = self.set_compound_het_var(var1, "111")
        var2 = self.set_compound_het_var(var2, "111")
        self.assertFalse(self.inh.is_compound_pair(var1, var2))
        
        
