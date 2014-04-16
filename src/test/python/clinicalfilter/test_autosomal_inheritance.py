""" unit testing of the Autosomal class
"""

import unittest
from clinicalfilter.ped import Family
from clinicalfilter.ped import Person
from clinicalfilter.variant import Variant
from clinicalfilter.variant_cnv import CNV
from clinicalfilter.variant_snv import SNV
from clinicalfilter.inheritance import Autosomal
from clinicalfilter.vcf_info import VcfInfo
from clinicalfilter.trio_genotypes import TrioGenotypes


class TestAutosomalPy(unittest.TestCase):
    """ test the Autosomal class
    """
    
    def setUp(self):
        """ define a family and variant, and start the Autosomal class
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
        
        var = TrioGenotypes(child_var)
        var.add_mother_variant(mom_var)
        var.add_father_variant(dad_var)
        self.variants = [var]
        
        # make sure we've got known genes data
        self.known_genes = {"TEST": {"inheritance": ["Monoallelic"], "confirmed_status": ["Confirmed DD Gene"]}}
        gene_inh = self.known_genes[var.get_gene()]["inheritance"]
        
        self.inh = Autosomal(self.variants, self.trio, gene_inh)
        self.inh.is_lof = var.child.is_lof()
    
    def create_snv(self, gender, genotype):
        """ create a default variant
        """
        
        chrom = "1"
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "G"
        qual = "50"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = SNV(chrom, pos, snp_id, ref, alt, qual, filt)
        
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
        qual = "50"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = CNV(chrom, pos, snp_id, ref, alt, qual, filt)
        
        tags = {"gene": ["HGNC", "VGN", "GN"], "consequence": ["VCQ", "CQ"]}
        
        info = "HGNC=TEST;HGNC_ALL=TEST;END=16000000;SVLEN=5000"
        format_keys = "INHERITANCE:DP"
        sample_values = inh + ":50"
        
        var.add_info(info, tags)
        var.add_format(format_keys, sample_values)
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
        self.assertEqual(self.inh.log_string, "de novo")
        
        # check for biallelic inheritance
        self.assertEqual(self.inh.check_heterozygous("Biallelic"), "compound_het")
        self.assertEqual(self.inh.log_string, "de novo")
        
        
        for geno in ["101", "102", "110", "112", "122"]:
            self.set_trio_genos(var, geno)
            self.inh.check_heterozygous("Monoallelic")
            self.assertNotEqual(self.inh.log_string, "de novo")
        
    def test_check_heterozygous_affected_mother(self):
        """ test that check_heterozygous() works correctly for affected mothers
        """
        
        var = self.variants[0]
        
        # check that trio = 110, with het affected mother is captured
        self.set_trio_genos(var, "110")
        self.inh.mother_affected = True
        self.assertEqual(self.inh.check_heterozygous("Monoallelic"), "single_variant")
        self.assertEqual(self.inh.log_string, "transmitted from aff, other parent non-carrier or aff")
        
        # check that when the other parent is also non-ref, the variant is no 
        # longer captured, unless the parent is affected
        self.set_trio_genos(var, "111")
        self.assertEqual(self.inh.check_heterozygous("Monoallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "typically for trios with non-de novo unaffected parents")
        
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
        child_var = self.create_cnv("F", "unknown", chrom, position)
        mom_var = self.create_cnv("F", "unknown", chrom, position)
        dad_var = self.create_cnv("M", "unknown", chrom, position)
        
        cnv_var = TrioGenotypes(child_var)
        cnv_var.add_mother_variant(mom_var)
        cnv_var.add_father_variant(dad_var)
        
        var = self.variants[0]
        
        # check for trio = 200, which is non-mendelian
        self.set_trio_genos(var, "200")
        self.assertEqual(self.inh.check_homozygous("Biallelic"), "nothing")
        self.assertEqual(self.inh.log_string, "non-mendelian trio")
        
        # check when a CNV is in the variants list
        self.inh.variants.append(cnv_var)
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


