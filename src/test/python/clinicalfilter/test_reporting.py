""" unit testing of the Report class
"""


import unittest
import logging
import os

from clinicalfilter.ped import Family
from clinicalfilter.ped import Person
from clinicalfilter.variant import Variant
from clinicalfilter.variant_cnv import CNV
from clinicalfilter.variant_snv import SNV
from clinicalfilter.trio_genotypes import TrioGenotypes
from clinicalfilter.reporting import Report

logging.disable(logging.CRITICAL)

class TestReportPy(unittest.TestCase):
    """ test the Report class
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
        child_var = self.create_snv(child_gender, "0/1")
        mom_var = self.create_snv("F", "0/0")
        dad_var = self.create_snv("M", "0/0")
        
        var = TrioGenotypes(child_var)
        var.add_mother_variant(mom_var)
        var.add_father_variant(dad_var)
        self.variants = [var]
        
        self.report = Report(None, None, None, None)
    
    def create_snv(self, gender, genotype):
        """ create a default variant
        """
        
        chrom = "X"
        pos = "15000000"
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
    
    def create_family(self, child_gender, mom_aff, dad_aff):
        """ create a default family, with optional gender and parental statuses
        """
        
        fam = Family("test")
        fam.add_child("child", "child_vcf", "2", child_gender)
        fam.add_mother("mother", "mother_vcf", mom_aff, "2")
        fam.add_father("father", "father_vcf", dad_aff, "1")
        fam.set_child()
        
        return fam
    
    def test__get_provenance(self):
        """ check that _get_provenance() works correctly
        """
        
        prov = ["checksum_string", "sample.calls.date.vcf.gz", "2014-01-01"]
        member = "proband"
        
        self.assertEqual(self.report._get_provenance(prov, member), \
            ["##UberVCF_proband_Id=sample\n", \
            "##UberVCF_proband_Checksum=checksum_string\n", \
            "##UberVCF_proband_Basename=sample.calls.date.vcf.gz\n", \
            "##UberVCF_proband_Date=2014-01-01\n"])
    
    def test__get_vcf_export_path(self):
        """ check that _get_vcf_export_path() works correctly
        """
        
        self.report.export_vcf
    
    