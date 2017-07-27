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
from clinicalfilter.inheritance import Autosomal
from clinicalfilter.trio_genotypes import TrioGenotypes

from tests.utils import create_snv, create_cnv

class TestAutosomalPy(unittest.TestCase):
    """ test the Autosomal class
    """
    
    def setUp(self):
        """ define a family and variant, and start the Autosomal class
        """
        
        # generate a test family
        sex = "F"
        mom_aff = "1"
        dad_aff = "1"
        
        self.trio = self.create_family(sex, mom_aff, dad_aff)
        
        # generate a test variant
        child = create_snv(sex, "0/1")
        mom = create_snv("F", "0/0")
        dad = create_snv("M", "0/0")
        
        var = TrioGenotypes(child.get_chrom(), child.get_position(),child, mom, dad)
        self.variants = [var]
        
        # make sure we've got known genes data
        self.known_gene = {"inh": ["Monoallelic"], "confirmed_status": ["confirmed dd gene"]}
        
        self.inh = Autosomal(self.variants, self.trio, self.known_gene, "1001")
        self.inh.is_lof = var.child.is_lof()
    
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
    
    def test_check_variant_without_parents(self):
        """ test that check_variant_without_parents() works correctly
        """
        
        var = self.variants[0]
        self.set_trio_genos(var, "100")
        
        # remove the parents, so it appears the var lacks parental information
        self.inh.trio.mother = None
        self.inh.trio.father = None
        
        # check for monoallelic inheritance
        self.assertEqual(self.inh.check_variant_without_parents("Monoallelic"), "single_variant")
        self.assertEqual(self.inh.log_string, "autosomal without parents")
        
        # check for biallelic inheritance
        self.assertEqual(self.inh.check_variant_without_parents("Biallelic"), "compound_het")
        
        self.set_trio_genos(var, "200")
        self.assertEqual(self.inh.check_variant_without_parents("Monoallelic"), "nothing")
        self.assertEqual(self.inh.check_variant_without_parents("Biallelic"), "single_variant")
    
    def test_check_heterozygous_de_novo(self):
        """ test that check_heterozygous() works correctly for de novos
        """
        
        var = self.variants[0]
        self.inh.set_trio_genotypes(var)
        # check for monoallelic inheritance
        self.assertEqual(self.inh.check_heterozygous("Monoallelic"), "single_variant")
        self.assertEqual(self.inh.log_string, "de novo as single_variant")
        
        # mosaic should operate as per monoallelic
        self.assertEqual(self.inh.check_heterozygous("Mosaic"), "single_variant")
        
        # check for biallelic inheritance
        self.assertEqual(self.inh.check_heterozygous("Biallelic"), "compound_het")
        self.assertEqual(self.inh.log_string, "de novo as compound_het")
        
        
        for geno in ["101", "102", "110", "112", "122"]:
            self.set_trio_genos(var, geno)
            self.inh.check_heterozygous("Monoallelic")
            self.assertNotEqual(self.inh.log_string, "de novo as single_variant")
        
    def test_check_heterozygous_affected_mother(self):
        """ test that check_heterozygous() works correctly for affected mothers
        """
        
        var = self.variants[0]
        
        # check that trio = 110, with het affected mother is captured
        self.set_trio_genos(var, "110")
        self.inh.mother_affected = True
        self.assertEqual(self.inh.check_heterozygous("Monoallelic"), "single_variant")
        self.assertEqual(self.inh.log_string, "transmitted from aff, other parent non-carrier or aff")
        
        # mosaic should operate as per monoallelic
        self.assertEqual(self.inh.check_heterozygous("Mosaic"), "single_variant")
        
        # check that when the other parent is also non-ref, the variant is no
        # longer captured, unless the parent is affected
        self.set_trio_genos(var, "111")
        self.assertEqual(self.inh.check_heterozygous("Monoallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "typically for trios with non-de novo unaffected parents")
        
        # mosaic should operate as per monoallelic
        self.assertEqual(self.inh.check_heterozygous("Mosaic"), "nothing")
        
        self.inh.father_affected = True
        self.inh.check_heterozygous("Monoallelic")
        self.assertEqual(self.inh.log_string, "transmitted from aff, other parent non-carrier or aff")
    
    def test_check_heterozygous_affected_father(self):
        """ test that check_heterozygous() works correctly for affected fathers
        """
        
        var = self.variants[0]
        # set the father as non-ref genotype and affected
        self.set_trio_genos(var, "101")
        self.inh.father_affected = True
        
        # check that the het proband, with het affected father is captured
        self.assertEqual(self.inh.check_heterozygous("Monoallelic"), "single_variant")
        self.assertEqual(self.inh.log_string, "transmitted from aff, other parent non-carrier or aff")
        
        # check that when the other parent is also non-ref, the variant is no
        # longer captured, unless the parent is affected
        self.set_trio_genos(var, "111")
        self.assertEqual(self.inh.check_heterozygous("Monoallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "typically for trios with non-de novo unaffected parents")
        
        self.inh.mother_affected = True
        self.inh.check_heterozygous("Monoallelic")
        self.assertEqual(self.inh.log_string, "transmitted from aff, other parent non-carrier or aff")
    
    def test_check_heterozygous_biallelic_non_ref_father(self):
        """ test that check_heterozygous() works correctly for biallelic genes
        """
        
        var = self.variants[0]
        # check for trio = 101
        self.set_trio_genos(var, "101")
        self.assertEqual(self.inh.check_heterozygous("Biallelic"), "compound_het")
        self.assertEqual(self.inh.log_string, "het-check for recessive genes and unaff parents not homoz")
        
        # check that trio = 102, with affected father comes through other route
        self.set_trio_genos(var, "102")
        self.inh.father_affected = True
        self.assertEqual(self.inh.check_heterozygous("Biallelic"), "compound_het")
        self.assertEqual(self.inh.log_string, "transmitted from aff, other parent non-carrier or aff")
        
        # check for trio = 102, without affected father
        self.inh.father_affected = False
        self.assertEqual(self.inh.check_heterozygous("Biallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "typically for trios with non-de novo unaffected parents")
    
    def test_check_heterozygous_biallelic_non_ref_mother(self):
        """ test that check_heterozygous() works correctly for biallelic genes
        """
        
        # check for trio = 110
        var = self.variants[0]
        self.set_trio_genos(var, "110")
        self.assertEqual(self.inh.check_heterozygous("Biallelic"), "compound_het")
        self.assertEqual(self.inh.log_string, "het-check for recessive genes and unaff parents not homoz")
        
        # check for trio = 120, with affected mother comes through other route
        var = self.variants[0]
        self.set_trio_genos(var, "120")
        self.inh.mother_affected = True
        self.assertEqual(self.inh.check_heterozygous("Biallelic"), "compound_het")
        self.assertEqual(self.inh.log_string, "transmitted from aff, other parent non-carrier or aff")
        
        # check for trio = 120, without affected mother
        self.set_trio_genos(var, "120")
        self.inh.mother_affected = False
        self.assertEqual(self.inh.check_heterozygous("Biallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "typically for trios with non-de novo unaffected parents")
    
    def test_check_heterozygous_biallelic_both_parents_non_ref(self):
        """ test that check_heterozygous() works correctly for biallelic genes
        """
        
        # check for trio = 111
        var = self.variants[0]
        self.set_trio_genos(var, "111")
        self.assertEqual(self.inh.check_heterozygous("Biallelic"), "compound_het")
        self.assertEqual(self.inh.log_string, "het-check for recessive genes and unaff parents not homoz")
        
        # check for trio = 121 with affected mother, but unaffected father
        self.set_trio_genos(var, "121")
        self.inh.mother_affected = True
        self.assertEqual(self.inh.check_heterozygous("Biallelic"), "compound_het")
        self.assertEqual(self.inh.log_string, "het-check for recessive genes and unaff parents not homoz")
        
        # check for trio = 112 with affected father, but unaffected mother
        self.set_trio_genos(var, "112")
        self.inh.mother_affected = False
        self.inh.father_affected = True
        
        self.assertEqual(self.inh.check_heterozygous("Biallelic"), "compound_het")
        self.assertEqual(self.inh.log_string, "het-check for recessive genes and unaff parents not homoz")
    
    def test_check_homozygous(self):
        """ test that check_homozygous() works correctly
        """
        
        var = self.variants[0]
        
        # check for trio = 200, which is non-mendelian
        self.set_trio_genos(var, "200")
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "non-mendelian trio")
        
        # check for 201, which is non-mendelian
        self.set_trio_genos(var, "201")
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "non-mendelian trio")
        
        # check for 210, which is non-mendelian
        self.set_trio_genos(var, "210")
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "non-mendelian trio")
        
    def test_check_homozygous_with_cnv(self):
        """ test that check_homozygous() works correctly for variant lists with CNVs
        """
        
        # generate a test variant
        chrom = "1"
        position = "60000"
        child = create_cnv("F", "unknown", chrom=chrom, pos=position)
        mom = create_cnv("F", "unknown", chrom=chrom, pos=position)
        dad = create_cnv("M", "unknown", chrom=chrom, pos=position)
        
        cnv = TrioGenotypes(chrom, position, child, mom, dad)
        var = self.variants[0]
        
        # check for trio = 200, which is non-mendelian
        self.set_trio_genos(var, "200")
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "non-mendelian trio")
        
        # check when a CNV is in the variants list
        self.inh.variants.append(cnv)
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "compound_het")
        self.assertEqual(self.inh.log_string, "non-mendelian, but CNV might affect call")
    
    def test_check_homozygous_biallelic(self):
        """ test that check_homozygous() works correctly for biallelic genes
        """
        
        var = self.variants[0]
        
        # check for trio = 211
        self.set_trio_genos(var, "211")
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "single_variant")
        self.assertEqual(self.inh.log_string, "both parents het in biallelic gene")
        
        # and check for trio = 221, but unaffected mother
        self.set_trio_genos(var, "221")
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "non-causal homozygous variant")
        
        # and check for trio = 221, with affected mother
        self.inh.mother_affected = True
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "single_variant")
        self.assertEqual(self.inh.log_string, "homoz parent aff")
        
        # and check for trio = 222, with affected mother only, should fail
        self.set_trio_genos(var, "222")
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "non-causal homozygous variant")
        
        # and check for trio = 222, with both parents affected
        self.inh.father_affected = True
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "single_variant")
        self.assertEqual(self.inh.log_string, "homoz parent aff")
        
        # and check for trio = 212, with affected father
        self.set_trio_genos(var, "212")
        self.inh.father_affected = True
        self.inh.mother_affected = False
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "single_variant")
        self.assertEqual(self.inh.log_string, "homoz parent aff")
    
    def test_check_homozygous_monoallelic(self):
        """ test that check_homozygous() works correctly for monoallelic genes
        """
        
        var = self.variants[0]
        
        # check for trio = 211, without affected father
        self.set_trio_genos(var, "211")
        self.assertEqual(self.inh.check_homozygous("Monoallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "non-causal homozygous variant")
        
        # check for trio = 211, with affected father only
        self.inh.father_affected = True
        self.assertEqual(self.inh.check_homozygous("Monoallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "non-causal homozygous variant")
        
        # check for trio = 211, with both parents affected
        self.inh.mother_affected = True
        self.assertEqual(self.inh.check_homozygous("Monoallelic"), "single_variant")
        self.assertEqual(self.inh.log_string, "transmitted from affected parents")
    


if __name__ == '__main__':
    unittest.main()
