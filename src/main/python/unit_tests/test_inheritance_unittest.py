""" unit testing of the Inheritance class
"""

import unittest
from clinicalfilter.ped import Family
from clinicalfilter.ped import Person
from clinicalfilter.variant import Variant
from clinicalfilter.variant_cnv import CNV
from clinicalfilter.variant_snv import SNV
from clinicalfilter.inheritance import Inheritance
from clinicalfilter.vcf_info import VcfInfo
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
        
        self.inh = Inheritance(self.variants, self.trio, gene_inh)
    
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
        self.format_keys = "GT:DP"
        self.sample_values = genotype + ":50"
        
        var.add_info(info, tags)
        var.set_gender(gender)
        
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
        
        self.assertTrue(self.inh.check_inheritance_mode_matches_gene_mode())
        
        


unittest.main()


        