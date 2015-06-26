""" unit testing of the VcfInfo class
"""

import unittest
from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.trio_genotypes import TrioGenotypes
from clinicalfilter.post_inheritance_filter import PostInheritanceFilter

class TestPostInheritanceFilterPy(unittest.TestCase):
    """
    """
    
    def setUp(self):
        """ define a default VcfInfo object
        """
        
        variants = []
        snv = self.create_var("1", True)
        cnv = self.create_var("1", False)
        
        variants.append((snv, "single_variant", "Monoallelic", ["ATRX"]))
        variants.append((cnv, "single_variant", "Monoallelic", ["ATRX"]))
        
        self.post_filter = PostInheritanceFilter(variants)
        
    def create_var(self, chrom, snv=True, geno=["0/1", "0/1", "0/1"]):
        """ define a family and variant, and start the Inheritance class
        
        Args:
            chrom: string for chrom, since we check the number of different chroms
            snv: boolean for whether to create a SNV or CNV object
        """
        
        # generate a test variant
        if snv:
            child_var = self.create_snv(chrom, geno[0])
            mom_var = self.create_snv(chrom, geno[1])
            dad_var = self.create_snv(chrom, geno[2])
        else:
            child_var = self.create_cnv(chrom)
            mom_var = self.create_cnv(chrom)
            dad_var = self.create_cnv(chrom)
        
        var = TrioGenotypes(child_var)
        var.add_mother_variant(mom_var)
        var.add_father_variant(dad_var)
        
        return var
    
    def create_snv(self, chrom, geno="0/1"):
        """ create a default variant
        """
        
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "G"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = SNV(chrom, pos, snp_id, ref, alt, filt)
        
        default_info = "HGNC=ATRX;CQ=missense_variant;random_tag;AF_AFR=0.0001"
        keys = "GT:DP:TEAM29_FILTER:PP_DNM"
        values = "{0}:50:PASS:0.99".format(geno)
        
        var.add_info(default_info)
        var.add_format(keys, values)
        var.set_gender("male")
        var.set_genotype()
        
        return var
    
    def create_cnv(self, chrom):
        
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "<DUP>"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = CNV(chrom, pos, snp_id, ref, alt, filt)
        
        info = "HGNC=TEST;HGNC_ALL=TEST,OR5A1;CQ=missense_variant;CNSOLIDATE;WSCORE=0.5;CALLP=0.000;COMMONFORWARDS=0.000;MEANLR2=0.5;MADL2R=0.02;END=16000000;SVLEN=1000000"
        format_keys = "inheritance:DP"
        sample_values = "deNovo:50"
        
        var.add_info(info)
        var.add_format(format_keys, sample_values)
        var.set_gender("F")
        var.set_genotype()
        
        return var
    
    def test_filter_variants(self):
        """ test that filter_variants() works correctly
        
        We only need to make a quick check of the first couple of lines, since
        the called functions are themselves tested elsewhere
        """
        
        variants = [(self.create_var("1", snv=False), "single_variant", "Biallelic", ["ATRX"])]
        variants.append((self.create_var("2", snv=False), "single_variant", "Biallelic", ["ATRX"]))
        
        # check that if we have CNVs on two chroms pass the filter
        # self.post_filter.variants = variants
        # self.assertEqual(self.post_filter.filter_variants(), variants)
        
        # check that CNVs on three different chroms get filtered out
        variants.append((self.create_var("3", snv=False), "single_variant", "Biallelic", ["ATRX"]))
        self.post_filter.variants = variants
        self.assertEqual(self.post_filter.filter_variants(), [])
    
    def test_count_cnv_chroms(self):
        """ test that count_cnv_chroms() works correctly
        """
        
        # check that the default list of variants counts only one chrom
        variants = self.post_filter.variants
        self.assertEqual(self.post_filter.count_cnv_chroms(variants), 1)
        
        # add CNVs on the same chrom, and check the chrom count increments one
        chrom_2_cnv_1 = self.create_var("2", snv=False)
        chrom_2_cnv_2 = self.create_var("2", snv=False)
        chrom_2_cnv_3 = self.create_var("2", snv=False)
        variants.append((chrom_2_cnv_1, "single_variant", "Biallelic", ["ATRX"]))
        variants.append((chrom_2_cnv_2, "single_variant", "Biallelic", ["ATRX"]))
        variants.append((chrom_2_cnv_3, "single_variant", "Biallelic", ["ATRX"]))
        self.assertEqual(self.post_filter.count_cnv_chroms(variants), 2)
        
        # and a CNV on a third chrom makes three
        chrom_3_cnv = self.create_var("3", snv=False)
        variants.append((chrom_3_cnv, "single_variant", "Biallelic", ["ATRX"]))
        self.assertEqual(self.post_filter.count_cnv_chroms(variants), 3)
    
    def test_remove_cnvs(self):
        """ test that remove_cnvs() works correctly
        """
        
        mixed_list = self.post_filter.variants
        snv_list = [mixed_list[0]]
        cnv_list = [mixed_list[1]]
        
        # check the effect of removing CNVs on different combinations of vars
        self.assertEqual(self.post_filter.remove_cnvs(mixed_list), snv_list)
        self.assertEqual(self.post_filter.remove_cnvs(snv_list), snv_list)
        self.assertEqual(self.post_filter.remove_cnvs(cnv_list), [])
    
    def test_filter_by_maf(self):
        """ test that filter_by_maf() works correctly
        """
        
        snv_1 = self.create_var("1", snv=True)
        snv_2 = self.create_var("2", snv=True)
        
        snv_1.child.info["AFR_AF"] = 0.0001
        snv_2.child.info["AFR_AF"] = 0.002
        
        # low maf Biallelic var returns the same
        variants = [(snv_1, "single_variant", "Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_by_maf(variants), variants)
        
        # low maf non-biallelic var returns the same
        variants = [(snv_1, "single_variant", "Monoallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_by_maf(variants), variants)
        
        # high maf Biallelic var returns the same
        variants = [(snv_2, "single_variant", "Monoallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_by_maf(variants), [])
        
        # high maf non-Biallelic is filtered out
        variants = [(snv_2, "single_variant", "Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_by_maf(variants), variants)
        
        # var with multiple inheritance modes should drop the non-biallelic
        # mode if the var has a high maf (leaving the Biallelic mode)
        variants = [(snv_2, "single_variant", "Monoallelic,Biallelic", ["ATRX"])]
        expected = [(snv_2, "single_variant", "Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_by_maf(variants), expected)
        
        # var with multiple inheritance modes should keep the non-biallelic
        # mode if the var has a low maf
        variants = [(snv_1, "single_variant", "Monoallelic,Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_by_maf(variants), variants)
        
        # check a de novo (lacking any MAF values)
        del snv_1.child.info["AFR_AF"]
        variants = [(snv_1, "single_variant", "Monoallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_by_maf(variants), variants)
        
    def test_filter_polyphen(self):
        """ check that filter_polyphen() works correctly
        """
        
        snv_1 = self.create_var("1", snv=True, geno=["0/1", "0/0", "0/1"])
        snv_2 = self.create_var("1", snv=True, geno=["0/1", "1/0", "0/1"])
        snv_3 = self.create_var("1", snv=True, geno=["0/1", "0/0", "0/0"])
        snv_1.position = 1000
        snv_2.position = 2000
        snv_3.position = 3000
        
        variants = [(snv_1, "single_variant", "Biallelic", ["ATRX"]), \
            (snv_2, "single_variant", "Biallelic", ["ATRX"])]
        
        # check that two vars without polyphen predictions pass
        self.assertEqual(self.post_filter.filter_polyphen(variants), variants)
        
        # check that two compound_hets in the same gene, with polyphen benign,
        # fail to pass the filter
        snv_1.child.info["PolyPhen"] = "benign(0.01)"
        snv_2.child.info["PolyPhen"] = "benign(0.01)"
        variants = [(snv_1, "compound_het", "Biallelic", ["ATRX"]), \
            (snv_2, "compound_het", "Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_polyphen(variants), [])
        
        # check that if one var is not benign, both compound hets fail to pass
        # the filter
        snv_2.child.info["PolyPhen"] = "probably_damaging(0.99)"
        self.assertEqual(self.post_filter.filter_polyphen(variants), [])
        
        # check that if one var lacks a polyphen value, and the other is damaging
        # both vars pass
        del snv_1.child.info["PolyPhen"]
        self.assertEqual(self.post_filter.filter_polyphen(variants), variants)
        
        # check that if both vars lack polyphen values, both vars pass
        del snv_2.child.info["PolyPhen"]
        self.assertEqual(self.post_filter.filter_polyphen(variants), variants)
        
        # check that single vars with polyphen benign fail
        snv_1.child.info["PolyPhen"] = "benign"
        variants = [(snv_1, "single_variant", "Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_polyphen(variants), [])
        
        # check if we have three compound_hets in the same gene, and one is
        # polyphen not benign, then all compound hets in the gene still pass,
        # even if two of them have polyphen benign
        snv_2.child.info["PolyPhen"] = "benign(0.01)"
        snv_3.child.info["PolyPhen"] = "probably_damaging(0.99)"
        variants = [(snv_1, "compound_het", "Biallelic", ["ATRX"]), \
            (snv_2, "compound_het", "Biallelic", ["ATRX"]), \
            (snv_3, "compound_het", "Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_polyphen(variants), [])
        
        # check if we have three compound_hets in the same gene, and two are
        # polyphen not benign, then only the two not benign compound hets pass
        snv_2.child.info["PolyPhen"] = "probably_damaging(0.99)"
        passing_vars = [(snv_2, "compound_het", "Biallelic", ["ATRX"]), \
            (snv_3, "compound_het", "Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_polyphen(variants), passing_vars)
    
    def test_has_compound_match(self):
        """ check that has_compound_match() works correctly
        """
        
        snv_1 = self.create_var("1", snv=True, geno=["0/1", "0/0", "0/1"])
        snv_2 = self.create_var("1", snv=True, geno=["0/1", "1/0", "0/1"])
        snv_3 = self.create_var("1", snv=True, geno=["0/1", "0/0", "0/0"])
        snv_1.position = 1000
        snv_2.position = 2000
        snv_3.position = 3000
        
        variants = [(snv_1, "compound_het", "Biallelic", ["ATRX"]), \
            (snv_2, "compound_het", "Biallelic", ["ATRX"])]
        
        # check that two vars without polyphen annotations return false
        self.assertFalse(self.post_filter.has_compound_match(snv_1, "ATRX", variants))
        
        # check that two vars with polyphen benign return true
        snv_1.child.info["PolyPhen"] = "benign(0.01)"
        snv_2.child.info["PolyPhen"] = "benign(0.01)"
        self.assertTrue(self.post_filter.has_compound_match(snv_1, "ATRX", variants))
        
        # check that having one var not polyphen benign returns True
        snv_2.child.info["PolyPhen"] = "probably_damaging(0.99)"
        self.assertTrue(self.post_filter.has_compound_match(snv_1, "ATRX", variants))
        
        # check that, if there are more than two compound hets to check in the
        # gene, we need two passing variants in order to pass
        snv_2.child.info["PolyPhen"] = "benign(0.01)"
        variants = [(snv_1, "compound_het", "Biallelic", ["ATRX"]), \
            (snv_2, "compound_het", "Biallelic", ["ATRX"]),
            (snv_3, "compound_het", "Biallelic", ["ATRX"])]
        self.assertTrue(self.post_filter.has_compound_match(snv_1, "ATRX", variants))
        
        # check that if we are checking a benign variant, and there are more
        # than two compound hets to check in the gene, if we have more than
        # two non-benign variants would prevent a match, then the function
        # returns false
        snv_1.child.info["PolyPhen"] = "probably_damaging(0.99)"
        variants = [(snv_1, "compound_het", "Biallelic", ["ATRX"]), \
            (snv_2, "compound_het", "Biallelic", ["ATRX"]),
            (snv_3, "compound_het", "Biallelic", ["ATRX"])]
        self.assertFalse(self.post_filter.has_compound_match(snv_1, "ATRX", variants))
        
        # check that we exclude benign de novos
        snv_1.child.info["PolyPhen"] = "probably_damaging(0.99)"
        snv_3.child.info["PolyPhen"] = "benign(0.01)"
        variants = [(snv_1, "compound_het", "Biallelic", ["ATRX"]), \
            (snv_3, "compound_het", "Biallelic", ["ATRX"])]
        self.assertFalse(self.post_filter.has_compound_match(snv_1, "ATRX", variants))
        
        # check that single variants in the same gene still return True
        variants = [(snv_1, "compound_het", "Biallelic", ["ATRX"]), \
            (snv_2, "single_variant", "Biallelic", ["ATRX"])]
        self.assertTrue(self.post_filter.has_compound_match(snv_1, "ATRX", variants))
    
    def test_filter_exac_hemizygous(self):
        """ check that filter_exac_hemizygous() works correctly
        """
        
        # construct a variant that will pass
        var = self.create_var("1", snv=True, geno=["0/1", "0/1", "0/1"])
        variants = [(var, "single_variant", "Biallelic", ["ATRX"])]
        
        # we should get back the same list of variants, if none of them have a
        # male chrX
        self.assertEqual(self.post_filter.filter_exac_hemizygous(variants), variants)
        
        # now construct a male chrX variant, which contains a non-zero AC_Hemi
        # annotation. This should fail the filter
        var = self.create_var("X", snv=True, geno=["1/1", "1/1", "1/1"])
        var.child.info["AC_Hemi"] = 1
        variants = [(var, "single_variant", "Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_exac_hemizygous(variants), [])
        
        # if the AC_Hemi count is zero, this should pass the filter
        var.child.info["AC_Hemi"] = 0
        variants = [(var, "single_variant", "Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_exac_hemizygous(variants), variants)
        
        # checkthat chrX females with non-zero AC_Hemi counts are not excluded
        var.inheritance_type = "XChrFemale"
        var.child.info["AC_Hemi"] = 1
        variants = [(var, "single_variant", "Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_exac_hemizygous(variants), variants)
        
        # now construct a de novo male chrX variant, which contains a non-zero
        # AC_Hemi annotation. Since this is not inherited, it should pass.
        var = self.create_var("X", snv=True, geno=["1/1", "0/0", "0/0"])
        var.child.info["AC_Hemi"] = 1
        variants = [(var, "single_variant", "Biallelic", ["ATRX"])]
        self.assertEqual(self.post_filter.filter_exac_hemizygous(variants), variants)
