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
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.variant.snv import SNV
from clinicalfilter.inheritance import Allosomal
from clinicalfilter.trio_genotypes import TrioGenotypes

from tests.utils import create_snv, create_cnv

class TestAllosomalPy(unittest.TestCase):
    """ test the Allosomal class
    """
    
    def setUp(self):
        """ define a family and variant, and start the Allosomal class
        """
        
        # generate a test family
        child_gender = "F"
        mom_aff = "1"
        dad_aff = "1"
        
        self.trio = self.create_family(child_gender, mom_aff, dad_aff)
        
        # generate a test variant
        child = create_snv(child_gender, "0/1")
        mom = create_snv("F", "0/0")
        dad = create_snv("M", "0/0")
        
        var = TrioGenotypes(child.get_chrom(), child.get_position(), child, mom, dad)
        self.variants = [var]
        
        # make sure we've got known genes data
        self.known_gene = {"inh": ["Monoallelic"], "confirmed_status": ["confirmed dd gene"]}
        
        self.inh = Allosomal(self.variants, self.trio, self.known_gene, "TEST")
        self.inh.is_lof = var.child.is_lof()
    
    def create_snv(self, gender, genotype):
        """ create a default variant
        """
        
        chrom = "X"
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "G"
        filt = "PASS"
        
        info = "HGNC=TEST;CQ=missense_variant;random_tag"
        format_keys = "GT:DP"
        sample_values = genotype + ":50"
        
        # set up a SNV object, since SNV inherits VcfInfo
        return SNV(chrom, pos, snp_id, ref, alt, filt, info, format_keys,
            sample_values, gender)
        
    def create_cnv(self, gender, inh, chrom, pos):
        """ create a default variant
        """
        
        snp_id = "."
        ref = "A"
        alt = "<DUP>"
        filt = "PASS"
        
        info = "HGNC=TEST;HGNC_ALL=TEST;END=16000000;SVLEN=5000"
        format_keys = "INHERITANCE:DP"
        sample_values = inh + ":50"
        
        # set up a SNV object, since SNV inherits VcfInfo
        return CNV(chrom, pos, snp_id, ref, alt, filt, info, format_keys,
            sample_values, gender)
        
    
    def create_family(self, child_gender, mom_aff, dad_aff):
        """ create a default family, with optional gender and parental statuses
        """
        
        fam = Family('test')
        fam.add_child('child', 'mother', 'father', child_gender, '2', 'child_vcf')
        fam.add_mother('mother', '0', '0', 'female', mom_aff, 'mother_vcf')
        fam.add_father('father', '0', '0', 'male', dad_aff, 'father_vcf')
        fam.set_child()
        
        return fam
    
    def set_trio_genos(self, var, geno):
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
        self.inh.set_trio_genotypes(var)
    
    def test_check_variant_without_parents_female(self):
        """ test that check_variant_without_parents() works correctly for female
        """
        
        var = TrioGenotypes('X', 100, self.create_snv('F', "1/0"), None, None)
        self.inh.set_trio_genotypes(var)
        
        # check for X-linked dominant inheritance
        self.assertEqual(self.inh.check_variant_without_parents("X-linked dominant"), "single_variant")
        self.assertEqual(self.inh.log_string, "allosomal without parents")
        
        var = TrioGenotypes('X', 100, self.create_snv('F', "1/1"),
            self.create_snv('F', "0/0"), self.create_snv('M', "0/0"))
        self.inh.set_trio_genotypes(var)
        self.assertEqual(self.inh.check_variant_without_parents("X-linked dominant"), "single_variant")
        
        # and check for hemizygous inheritance
        var = TrioGenotypes('X', 100, self.create_snv('F', "1/0"),
            self.create_snv('F', "0/0"), self.create_snv('M', "0/0"))
        self.inh.set_trio_genotypes(var)
        self.assertEqual(self.inh.check_variant_without_parents("Hemizygous"), "hemizygous")
        
        var = TrioGenotypes('X', 100, self.create_snv('F', "1/1"),
            self.create_snv('F', "0/0"), self.create_snv('M', "0/0"))
        self.inh.set_trio_genotypes(var)
        self.assertEqual(self.inh.check_variant_without_parents("Hemizygous"), "single_variant")
    
    def test_check_variant_without_parents_male(self):
        """ test that check_variant_without_parents() works correctly for males
        """
        
        var = TrioGenotypes('X', 100, self.create_snv('M', "1/1"), None, None)
        trio = self.create_family('male', '1', '1')
        
        self.inh = Allosomal([var], trio, self.known_gene, "TEST")
        self.inh.set_trio_genotypes(var)
        
        # check for X-linked dominant inheritance
        self.assertEqual(self.inh.check_variant_without_parents("X-linked dominant"), "single_variant")
        
        # and check for hemizygous inheritance
        self.assertEqual(self.inh.check_variant_without_parents("Hemizygous"), "single_variant")
    
    def test_check_heterozygous_de_novo(self):
        """ test that check_heterozygous() works correctly for de novos
        """
        
        # all of these tests are run for female X chrom de novos, since male
        # X chrom hets don't exist
        var = TrioGenotypes('X', 100, self.create_snv('F', "1/0"),
            self.create_snv('F', "0/0"), self.create_snv('M', "0/0"))
        self.inh.set_trio_genotypes(var)
        
        # check for X-linked dominant inheritance
        self.assertEqual(self.inh.check_heterozygous("X-linked dominant"), "single_variant")
        self.assertEqual(self.inh.log_string, "female x chrom de novo")
        
        # check for biallelic inheritance
        self.assertEqual(self.inh.check_heterozygous("Hemizygous"), "single_variant")
        
        # check for X-linked over dominance
        self.assertEqual(self.inh.check_heterozygous("X-linked over-dominance"), 'single_variant')
        
        # check we raise errors with unknown inheritance modes
        with self.assertRaises(ValueError):
            self.inh.check_heterozygous("Digenic")
        
        genos = {'0': '0/0', '1': '1/0', '2': '1/1'}
        
        for (child, mom, dad) in ["102", "110", "112", "122"]:
            var = TrioGenotypes('X', 100, self.create_snv('F', genos[child]),
                self.create_snv('F', genos[mom]), self.create_snv('M', genos[dad]))
            self.inh.set_trio_genotypes(var)
            
            self.inh.check_heterozygous("X-linked dominant")
            self.assertNotEqual(self.inh.log_string, "female x chrom de novo")
        
    def test_check_heterozygous_affected_mother(self):
        """ test that check_heterozygous() works correctly for affected mothers
        """
        
        # check that trio with het affected mother is captured
        var = TrioGenotypes('X', 100, self.create_snv('F', "1/0"),
            self.create_snv('F', "1/0"), self.create_snv('M', "0/0"))
        self.inh.set_trio_genotypes(var)
        
        self.inh.mother_affected = True
        self.assertEqual(self.inh.check_heterozygous("X-linked dominant"), "single_variant")
        self.assertEqual(self.inh.log_string, "x chrom transmitted from aff, other parent non-carrier or aff")
        
        # check that when the other parent is also non-ref, the variant is no
        # longer captured, unless the parent is affected
        var = TrioGenotypes('X', 100, self.create_snv('F', "1/0"),
            self.create_snv('F', "1/0"), self.create_snv('M', "1/1"))
        self.inh.set_trio_genotypes(var)
        self.assertEqual(self.inh.check_heterozygous("X-linked dominant"), "nothing")
        self.assertEqual(self.inh.log_string, "variant not compatible with being causal")
        
        self.inh.father_affected = True
        self.inh.check_heterozygous("X-linked dominant")
        self.assertEqual(self.inh.log_string, "x chrom transmitted from aff, other parent non-carrier or aff")
        
        # and check that hemizgygous vars return as "compound_het"
        self.assertEqual(self.inh.check_heterozygous("Hemizygous"), "compound_het")
    
    def test_check_heterozygous_affected_father(self):
        """ test that check_heterozygous() works correctly for affected fathers
        """
        
        # set the father as non-ref genotype and unaffected
        var = TrioGenotypes('X', 100, self.create_snv('F', "1/0"),
            self.create_snv('F', "0/0"), self.create_snv('M', "1/1"))
        self.inh.set_trio_genotypes(var)
        
        # check that the het proband, with het unaffected father is passes
        self.assertEqual(self.inh.check_heterozygous("X-linked dominant"), "nothing")
        self.assertEqual(self.inh.log_string, "variant not compatible with being causal")
        
        # check that the het proband, with het affected father is captured
        self.inh.father_affected = True
        self.assertEqual(self.inh.check_heterozygous("X-linked dominant"), "single_variant")
        self.assertEqual(self.inh.log_string, "x chrom transmitted from aff, other parent non-carrier or aff")
        
        # check that when the other parent is also non-ref, the variant is no
        # longer captured, unless the parent is affected
        var = TrioGenotypes('X', 100, self.create_snv('F', "1/0"),
            self.create_snv('F', "1/0"), self.create_snv('M', "1/1"))
        self.inh.set_trio_genotypes(var)
        
        self.assertEqual(self.inh.check_heterozygous("X-linked dominant"), "nothing")
        self.assertEqual(self.inh.log_string, "variant not compatible with being causal")
        
        self.inh.mother_affected = True
        self.inh.check_heterozygous("X-linked dominant")
        self.assertEqual(self.inh.log_string, "x chrom transmitted from aff, other parent non-carrier or aff")
    
    def test_check_homozygous_male(self):
        """ test that check_homozygous() works correctly for males
        """
        
        # check for trio with de novo on male X chrom
        var = TrioGenotypes('X', 100, self.create_snv('M', "1/1"),
            self.create_snv('F', "0/0"), self.create_snv('M', "0/0"))
        
        trio = self.create_family('male', '1', '1')
        self.inh = Allosomal([var], trio, self.known_gene, "TEST")
        self.inh.set_trio_genotypes(var)
        
        self.assertEqual(self.inh.check_homozygous("X-linked dominant"), "single_variant")
        self.assertEqual(self.inh.log_string, "male X chrom de novo")
        
        # check for trio = 210, with unaffected mother
        var = TrioGenotypes('X', 100, self.create_snv('M', "1/1"),
            self.create_snv('F', "1/0"), self.create_snv('M', "0/0"))
        self.inh.set_trio_genotypes(var)
        
        self.assertEqual(self.inh.check_homozygous("X-linked dominant"), "single_variant")
        self.assertEqual(self.inh.log_string, "male X chrom inherited from het mother or hom affected mother")
        
        # check for trio = 210, with affected mother, which should not pass
        self.inh.mother_affected = True
        self.assertEqual(self.inh.check_homozygous("X-linked dominant"), "nothing")
        self.assertEqual(self.inh.log_string, "variant not compatible with being causal")
        
        # check for trio = 220, with affected mother
        var = TrioGenotypes('X', 100, self.create_snv('M', "1/1"),
            self.create_snv('F', "1/1"), self.create_snv('M', "0/0"))
        self.inh.set_trio_genotypes(var)
        self.assertEqual(self.inh.check_homozygous("X-linked dominant"), "single_variant")
        self.assertEqual(self.inh.log_string, "male X chrom inherited from het mother or hom affected mother")
        
        # check for trio = 220, with unaffected mother, which should not pass
        self.inh.mother_affected = False
        self.assertEqual(self.inh.check_homozygous("X-linked dominant"), "nothing")
        self.assertEqual(self.inh.log_string, "variant not compatible with being causal")
        
        # check that homozygous X-linked over-dominance doesn't pass
        self.assertEqual(self.inh.check_homozygous("X-linked over-dominance"), 'nothing')
        
        # check we raise errors with unknown inheritance modes
        with self.assertRaises(ValueError):
            self.inh.check_homozygous("Digenic")
    
    def test_check_homozygous_female(self):
        """ test that check_homozygous() works correctly for females
        """
        
        var = self.variants[0]
        
        # check for trio = 200, which is non-mendelian
        self.set_trio_genos(var, "200")
        self.assertEqual(self.inh.check_homozygous("Hemizygous"), "nothing")
        self.assertEqual(self.inh.log_string, "non-mendelian trio")
        
        # check for trio = 210, which is non-mendelian
        self.set_trio_genos(var, "210")
        self.assertEqual(self.inh.check_homozygous("Hemizygous"), "nothing")
        self.assertEqual(self.inh.log_string, "non-mendelian trio")
        
        # check for trio = 202, which is non-mendelian
        self.set_trio_genos(var, "202")
        self.assertEqual(self.inh.check_homozygous("Hemizygous"), "nothing")
        self.assertEqual(self.inh.log_string, "non-mendelian trio")
        
        # and check for trio = 212, without affected parents
        self.set_trio_genos(var, "212")
        self.assertEqual(self.inh.check_homozygous("Hemizygous"), "nothing")
        self.assertEqual(self.inh.log_string, "variant not compatible with being causal")
        
        # and check for trio = 212, with affected father
        self.set_trio_genos(var, "212")
        self.inh.father_affected = True
        self.assertEqual(self.inh.check_homozygous("Hemizygous"), "single_variant")
        self.assertEqual(self.inh.log_string, "testing")
        
        # and check for trio = 212, with affected mother
        self.set_trio_genos(var, "212")
        self.inh.mother_affected = True
        self.assertEqual(self.inh.check_homozygous("Hemizygous"), "single_variant")
        self.assertEqual(self.inh.log_string, "testing")
        
        # and check for trio = 222, with affected mother
        self.set_trio_genos(var, "222")
        self.assertEqual(self.inh.check_homozygous("Hemizygous"), "single_variant")
        self.assertEqual(self.inh.log_string, "testing")
        
        # and check for trio = 222, with affected mother
        self.set_trio_genos(var, "222")
        self.inh.mother_affected = False
        self.assertEqual(self.inh.check_homozygous("Hemizygous"), "nothing")
        self.assertEqual(self.inh.log_string, "variant not compatible with being causal")
    
    def test_check_homozygous_with_cnv(self):
        """ test that check_homozygous() works correctly for variant lists with CNVs
        """
        
        # generate a test variant
        chrom = "X"
        pos = '160'
        child = create_cnv('F', 'unknown', chrom=chrom, pos=pos)
        mom = create_cnv('F', 'unknown', chrom=chrom, pos=pos)
        dad = create_cnv('F', 'unknown', chrom=chrom, pos=pos)
        
        cnv = TrioGenotypes(chrom, pos, child, mom, dad)
        
        var = self.variants[0]
        
        # check for trio = 200, which is non-mendelian
        self.set_trio_genos(var, "200")
        self.assertEqual(self.inh.check_homozygous("Hemizygous"), "nothing")
        self.assertEqual(self.inh.log_string, "non-mendelian trio")
        
        # check when a CNV is in the variants list
        self.inh.variants.append(cnv)
        self.assertEqual(self.inh.check_homozygous("Hemizygous"), "compound_het")
        self.assertEqual(self.inh.log_string, "non-mendelian, but CNV might affect call")
