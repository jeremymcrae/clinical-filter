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

from clinicalfilter.ped import Family
from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.variant.info import Info
from clinicalfilter.inheritance import Autosomal
from clinicalfilter.inheritance import Allosomal
from clinicalfilter.trio_genotypes import TrioGenotypes

from tests.utils import create_snv, create_cnv

class TestInheritancePy(unittest.TestCase):
    """ test the Inheritance class
    """
    
    def setUp(self):
        """ define a family and variant, and start the Inheritance class
        """
        
        # generate a test family
        sex = "F"
        mom_aff = "1"
        dad_aff = "1"
        
        self.trio = self.create_family(sex, mom_aff, dad_aff)
        
        # generate list of variants
        self.variants = [self.create_variant(sex)]
        self.variants.append(self.create_variant(sex))
        
        # make sure we've got known genes data
        self.known_gene = {"inh": ["Monoallelic"], "confirmed_status": ["confirmed dd gene"]}
        
        self.inh = Autosomal(self.variants, self.trio, self.known_gene, "1001")
    
    def create_variant(self, chrom="1", position="150", sex='F', cq=None,
            geno=['0/1', '0/0', '0/0']):
        """ creates a TrioGenotypes variant
        """
        
        # generate a test variant
        child = create_snv(sex, geno[0], cq, chrom=chrom, pos=position)
        mom = create_snv("F", geno[1], cq, chrom=chrom, pos=position)
        dad = create_snv("M", geno[2], cq, chrom=chrom, pos=position)
        
        return TrioGenotypes(chrom, position, child, mom, dad)
    
    def create_family(self, child_gender, mom_aff, dad_aff):
        """ create a default family, with optional gender and parental statuses
        """
        
        fam = Family('test')
        fam.add_child('child', 'mother', 'father', child_gender, '2', 'child_vcf')
        fam.add_mother('mother', '0', '0', 'female', mom_aff, 'mother_vcf')
        fam.add_father('father', '0', '0', 'male', dad_aff, 'father_vcf')
        fam.set_child()
        
        return fam
    
    def test_check_inheritance_mode_matches_gene_mode(self):
        """ test that check_inheritance_mode_matches_gene_mode() works correctly
        """
        
        # check that the default inheritance types have been set up correctly
        self.assertEqual(self.inh.inheritance_modes, {"Monoallelic",
            "Biallelic", "Both", 'Imprinted'})
        
        # make sure that the default var and gene inheritance work
        self.assertTrue(self.inh.check_inheritance_mode_matches_gene_mode())
        
        # check that no gene inheritance overlap fails
        self.inh.gene_inheritance = {"Mosaic"}
        self.inh.inheritance_modes = {"Monoallelic", "Biallelic", "Both", 'Imprinted'}
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
    
    def test_get_candidate_variants_monoallelic(self):
        """ test that get_candidate_variants() works for a monoallelic variant
        """
        
        inh = {"inh": ["Monoallelic"], "confirmed_status": ["confirmed dd gene"]}
        var = self.create_variant(position='150', cq='stop_gained',
            geno=['0/1', '0/0', '0/0'])
        self.inh = Autosomal([var], self.trio, inh, "TEST")
        
        self.assertEqual(self.inh.get_candidate_variants(),
            [(var, ['single_variant'], ['Monoallelic'], ['TEST'])])
        
        # check a variant that shouldn't pass the monoallelic route
        var = self.create_variant(position='150', cq='stop_gained',
            geno=['0/1', '0/1', '0/0'])
        self.inh = Autosomal([var], self.trio, inh, "TEST")
        self.assertEqual(self.inh.get_candidate_variants(), [])
    
    
    def test_get_candidate_variants_imprinted(self):
        """ test that get_candidate_variants() works for imprinted variants
        """
        
        # check a variant where the imprinting route should work
        inh = {"inh": ["Imprinted"], "confirmed_status": ["confirmed dd gene"]}
        var = self.create_variant(position='150', cq='stop_gained',
            geno=['0/1', '0/1', '0/0'])
        self.inh = Autosomal([var], self.trio, inh, "TEST")
        
        self.assertEqual(self.inh.get_candidate_variants(),
            [(var, ['single_variant'], ['Imprinted'], ['TEST'])])
        
        # de novos shouldn't pass the imprinted route
        var = self.create_variant(position='150', cq='stop_gained',
            geno=['0/1', '0/0', '0/0'])
        self.inh = Autosomal([var], self.trio, inh, "TEST")
        self.assertEqual(self.inh.get_candidate_variants(), [])
        
        # check a variant that shouldn't pass the imprinted route due to there
        # not being a known gene.
        # NOTE: this behavior differs differs slightly from other inheritance
        # modes when there isn't any known gene. Since there are so few known
        # imprinting genes, it makes sense to only allow for these when we know
        # the mode is correct
        var = self.create_variant(position='150', cq='stop_gained',
            geno=['0/1', '0/1', '0/0'])
        self.inh = Autosomal([var], self.trio, None, "TEST")
        self.assertEqual(self.inh.get_candidate_variants(), [])
        
        # also check imprinting requires loss-of-function consequences
        var = self.create_variant(position='150', cq='missense_variant',
            geno=['0/1', '0/1', '0/0'])
        self.inh = Autosomal([var], self.trio, inh, "TEST")
        self.assertEqual(self.inh.get_candidate_variants(), [])
        
        # check loss-of-function requirement for a paternally inherited variant
        var = self.create_variant(position='150', cq='missense_variant',
            geno=['0/1', '0/0', '0/1'])
        self.inh = Autosomal([var], self.trio, inh, "TEST")
        self.assertEqual(self.inh.get_candidate_variants(), [])
        
        # check imprinting for a biallelic variant
        var = self.create_variant(position='150', cq='stop_gained',
            geno=['1/1', '0/1', '0/1'])
        self.inh = Autosomal([var], self.trio, inh, "TEST")
        self.assertEqual(self.inh.get_candidate_variants(),
            [(var, ['single_variant'], ['Imprinted'], ['TEST'])])
        
        # check imprinting for a biallelic variant
        var = self.create_variant(position='150', cq='missense_variant',
            geno=['1/1', '0/1', '0/1'])
        self.inh = Autosomal([var], self.trio, inh, "TEST")
        self.assertEqual(self.inh.get_candidate_variants(), [])
    
    def test_get_candidate_variants_compound_het(self):
        """ test that get_candidate_variants() works for biallelic variants
        """
        
        inh = {"inh": ["Biallelic"], "confirmed_status": ["confirmed dd gene"]}
        var1 = self.create_variant(position='150', cq='stop_gained',
            geno=['0/1', '0/1', '0/0'])
        var2 = self.create_variant(position='151', cq='stop_gained',
            geno=['0/1', '0/0', '1/0'])
        self.inh = Autosomal([var1, var2], self.trio, inh, "TEST")
        
        self.assertEqual(sorted(self.inh.get_candidate_variants()),
            sorted([(var1, ['compound_het'], ['Biallelic'], ['TEST']),
            (var2, ['compound_het'], ['Biallelic'], ['TEST'])]))
        
        # check that a single variant isn't included in the compound hets
        self.inh = Autosomal([var1], self.trio, inh, "TEST")
        self.assertEqual(self.inh.get_candidate_variants(), [])
    
    def test_check_if_any_variant_is_cnv(self):
        """ test if check_if_any_variant_is_cnv() works correctly
        """
        
        # generate a test variant
        chrom = "1"
        position = "60000"
        child = create_cnv("F", "unknown", chrom=chrom, pos=position)
        mom = create_cnv("F", "unknown", chrom=chrom, pos=position)
        dad = create_cnv("M", "unknown", chrom=chrom, pos=position)
        
        cnv_var = TrioGenotypes(chrom, position, child, mom, dad)
        
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
        var1 = self.create_variant(chrom="1", position="150", sex="F", cq="stop_gained")
        var2 = self.create_variant(chrom="1", position="160", sex="F", cq="stop_gained")
        var3 = self.create_variant(chrom="1", position="170", sex="F", cq="stop_gained")
        
        # set the inheritance type, the compound het type ("compound_het"
        # for autosomal variants, and start autosomal inheritance)
        # known_genes = "Biallelic"
        known_gene = {"inh": ["Biallelic"], "confirmed_status": ["confirmed dd gene"]}
        self.inh = Autosomal([var1, var2, var3], self.trio, known_gene, "TEST")
        
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
        var1 = self.create_variant(chrom="1", position="150", sex="F", cq="stop_gained")
        var2 = self.create_variant(chrom="1", position="160", sex="F", cq="stop_gained")
        
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
        var1 = self.create_variant(chrom="1", position="150", sex="F", cq="missense_variant")
        var2 = self.create_variant(chrom="1", position="160", sex="F", cq="missense_variant")
        var3 = self.create_variant(chrom="1", position="160", sex="F", cq="inframe_deletion")
        var4 = self.create_variant(chrom="1", position="160", sex="F", cq="stop_gained")
        
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
        var1 = self.create_variant(chrom="1", position="150", sex="F", cq="missense_variant")
        var2 = self.create_variant(chrom="1", position="160", sex="F", cq="missense_variant")
        var3 = self.create_variant(chrom="1", position="160", sex="F", cq="inframe_deletion")
        var4 = self.create_variant(chrom="1", position="160", sex="F", cq="stop_gained")
        
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
        var1 = self.create_variant(chrom="1", position="150", sex="F", cq="stop_gained")
        var2 = self.create_variant(chrom="1", position="160", sex="F", cq="stop_gained")
        
        var1 = self.set_compound_het_var(var1, "110")
        var2 = self.set_compound_het_var(var2, "101")
        
        var1.child.info = Info('CQ=missense_variant')
        var2.child.info = Info('CQ=missense_variant')
        var1.child.info.set_genes_and_consequence('1', 100, ('G', ), [])
        var2.child.info.set_genes_and_consequence('1', 100, ('G', ), [])
        
        # exclude pairs where both members are not loss-of-function
        self.assertFalse(self.inh.is_compound_pair(var1, var2))
    
    def test_is_compound_pair_cnv_paternal(self):
        """ check that is_compound_pair() includes pairs with CNVs
        """
        
        # generate a test variant
        chrom = "1"
        position = "60000"
        extra = [('CIFER_INHERITANCE', 'paternal')]
        child = create_cnv("F", "paternal", chrom=chrom, pos=position, format=extra)
        mom = create_cnv("F", "unknown", chrom=chrom, pos=position)
        dad = create_cnv("M", "unknown", chrom=chrom, pos=position)
        
        cnv = TrioGenotypes(chrom, position, child, mom, dad)
        
        # set some variants, so we can alter them later
        snv = self.create_variant(chrom="1", position="150", sex="F", cq="stop_gained")
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
        extra = [('CIFER_INHERITANCE', 'maternal')]
        child = create_cnv("F", "maternal", chrom=chrom, pos=position, format=extra)
        mom = create_cnv("F", "unknown", chrom=chrom, pos=position)
        dad = create_cnv("M", "unknown", chrom=chrom, pos=position)
        
        cnv = TrioGenotypes(chrom, position, child, mom, dad)
        
        # set some variants, so we can alter them later
        snv = self.create_variant(chrom="1", position="150", sex="F", cq="stop_gained")
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
        fam.add_child("child", 'dad_id', 'mom_id', 'F', '2',  "child_vcf")
        fam.set_child()
        
        # set some variants, so we can alter them later
        var1 = self.create_variant(chrom="1", position="150", sex="F", cq="stop_gained")
        var2 = self.create_variant(chrom="1", position="160", sex="F", cq="stop_gained")
        
        inh = Autosomal([var1, var2], fam, self.known_gene, "TEST")
        
        # check that a proband-only passes, regardless of the parental genotypes
        self.assertTrue(inh.is_compound_pair(var1, var2))
    
    def test_is_compound_pair_allosomal(self):
        """ check that is_compound_pair() works when the father is affected
        """
        
        # set some variants, so we can alter them later
        var1 = self.create_variant(chrom="X", position="150", sex="F", cq="stop_gained")
        var2 = self.create_variant(chrom="X", position="160", sex="F", cq="stop_gained")
        
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
        var1 = self.create_variant(chrom="1", position="150", sex="M", cq="stop_gained")
        var2 = self.create_variant(chrom="1", position="160", sex="M", cq="stop_gained")
        
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
