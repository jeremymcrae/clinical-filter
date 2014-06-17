""" unit testing of the LoadVCFs class
"""

import unittest
import gzip
import zlib
import os
import sys
import io
import shutil
import tempfile
import random
import hashlib

from clinicalfilter.variant_snv import SNV
from clinicalfilter.variant_cnv import CNV
from clinicalfilter.trio_genotypes import TrioGenotypes
from clinicalfilter.match_cnvs import MatchCNVs
from clinicalfilter.vcf import LoadVCFs
from clinicalfilter.ped import Family, Person

IS_PYTHON2 = sys.version_info[0] == 2
IS_PYTHON3 = sys.version_info[0] == 3


class TestLoadVCFsPy(unittest.TestCase):
    """
    """
    
    def setUp(self):
        """ define a default LoadVCFs object
        """
        
        counter = 0
        total_trios = 1
        filters = {}
        tags_dict = {}
        
        self.vcf_loader = LoadVCFs(counter, total_trios, filters, tags_dict)
        
        # make a temp directory for the cache file
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """ remove the temp directory once a test completes
        """
        
        shutil.rmtree(self.temp_dir)
    
    def make_minimal_vcf(self):
        """ construct the bare minimum of lines for a VCF file
        """
        
        header = []
        header.append("##fileformat=VCFv4.1\n")
        header.append("##fileDate=2014-01-01\n")
        header.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_id\n")
        
        variants = []
        variants.append("1\t100\t.\tT\tA\t1000\tPASS\t.\tGT\t0/1\n")
        variants.append("1\t200\t.\tT\tA\t1000\tPASS\t.\tGT\t0/1\n")
        
        vcf = header + variants
        
        return vcf
    
    def write_temp_vcf(self, filename, vcf_data):
        """ writes data to a file, and returns the full path to the file
        """
        
        full_path = os.path.join(self.temp_dir, filename)
        
        vcf_data = "".join(vcf_data)
        output = open(full_path, "w")
        output.write(vcf_data)
        output.close()
        
        return full_path
    
    def write_gzipped_vcf(self, filename, vcf_data):
        """ writes data to a gzip file, and returns the full path to the file
        """
        
        full_path = os.path.join(self.temp_dir, filename)
        
        vcf_data = "".join(vcf_data)
        if IS_PYTHON2:
            f = gzip.open(full_path, 'wb')
        elif IS_PYTHON3:
            f = gzip.open(full_path, 'wt')
        f.write(vcf_data)
        f.close()
        
        return full_path
    
    def test_open_vcf_file(self):
        """ test obtaining a file handle for the VCF
        """
        
        vcf = self.make_minimal_vcf()
        path = self.write_temp_vcf("temp.vcf", vcf)
        
        # check that plain VCF files can be loaded
        handle = self.vcf_loader.open_vcf_file(path)
        self.assertEqual(type(handle), io.TextIOWrapper)
        handle.close()
        
        # check that gzipped vcf files are handled correctly
        path = self.write_gzipped_vcf("temp.vcf.gz", vcf)
        
        handle = self.vcf_loader.open_vcf_file(path)
        if IS_PYTHON2:
            self.assertEqual(type(handle), gzip.GzipFile)
        elif IS_PYTHON3:
            self.assertEqual(type(handle), io.TextIOWrapper)
        handle.close()
        
        # make sure files that don't exists raise an error
        path = os.path.join(self.temp_dir, "zzz.txt")
        with self.assertRaises(OSError):
            self.vcf_loader.open_vcf_file(path)
        
        # check that files with unknown extensions raise errors
        path = self.write_temp_vcf("temp.zzz", vcf)
        with self.assertRaises(OSError):
            self.vcf_loader.open_vcf_file(path)
    
    def test_get_vcf_header(self):
        """ test that get_vcf_header() works correctly
        """
        
        vcf = self.make_minimal_vcf()
        path = self.write_temp_vcf("temp.vcf", vcf)
        
        header = self.vcf_loader.get_vcf_header(path)
        
        # check that the header is returned correctly
        self.assertEqual(header, vcf[:4])
    
    def test_exclude_header(self):
        """ test that exclude_header() works correctly
        """
        
        vcf = self.make_minimal_vcf()
        
        # make sure we drop the header, and only the header from the file
        # check this by reading the file, and making sure the first line 
        # is the line we expect from the VCF
        path = self.write_temp_vcf("temp.vcf", vcf)
        handler = open(path, "r")
        self.vcf_loader.exclude_header(handler)
        self.assertEqual(handler.readline(), vcf[4])
        handler.close()
        
        # also check for gzipped VCF files.
        path = self.write_gzipped_vcf("temp.vcf.gz", vcf)
        if IS_PYTHON2:
            handler = gzip.open(path, "r")
        elif IS_PYTHON3:
            handler = gzip.open(path, "rt")
        self.vcf_loader.exclude_header(handler)
        self.assertEqual(handler.readline(), vcf[4])
        handler.close()
    
    def test_add_single_variant(self):
        """ test that add_single_variant() works correctly
        """
        
        # the sub-functions are all tested elsewhere, this test merely checks
        # that valid variants are added to the variants list, and invalid
        # variants are passed over without being added to the variants list
        
        # set up an autosomal variant
        line = ["1", "100", ".", "T", "G", "1000", "PASS", ".", "GT", "0/1"]
        gender = "M"
        variant = SNV(*line[:6])
        
        # check that the variant is added to the variant list
        variants = []
        self.vcf_loader.add_single_variant(variants, variant, gender, line)
        self.assertEqual(variants, [variant])
        
        # set up an X-chrom male het
        line = ["X", "100", ".", "T", "G", "1000", "PASS", ".", "GT", "0/1"]
        variant = SNV(*line[:6])
        
        # check that the X-chrom male het is not added to the variant list
        variants = []
        self.vcf_loader.add_single_variant(variants, variant, gender, line)
        self.assertEqual(variants, [])
    
    def test_get_vcf_provenance(self):
        """ test that get_vcf_provenance() works correctly
        """
        
        vcf = self.make_minimal_vcf()
        vcf_string = "".join(vcf)
        if IS_PYTHON3:
            vcf_string = vcf_string.encode("utf-8")
        ungzipped_hash = hashlib.sha1(vcf_string).hexdigest()
        header = vcf[:4]
        
        path = self.write_temp_vcf("temp.vcf", vcf)
        
        # check that the file defs return correctly
        (checksum, basename, date) = self.vcf_loader.get_vcf_provenance(path)
        
        self.assertEqual(checksum, ungzipped_hash)
        self.assertEqual(basename, "temp.vcf")
        self.assertEqual(date, "2014-01-01")
        
        # now write a gzip file, and check that we get the correct hash
        path = self.write_gzipped_vcf("test.vcf.gz", vcf)
        handle = open(path, "rb")
        gzipped_hash = hashlib.sha1(handle.read()).hexdigest()
        handle.close()
        
        (checksum, basename, date) = self.vcf_loader.get_vcf_provenance(path)
        self.assertEqual(checksum, gzipped_hash)
        
        # check that when a fileDate isn't available in the VCf, we can pick 
        # the date from the path
        vcf.pop(1)
        path = self.write_temp_vcf("temp.file_process.2014-02-20.vcf", vcf)
        (checksum, basename, date) = self.vcf_loader.get_vcf_provenance(path)
        self.assertEqual(date, "2014-02-20")
    
    def test_construct_variant(self):
        """ test that construct_variant() works correctly
        """
        
        # check that construct variant works for SNVs
        line = ["1", "100", ".", "T", "G", "1000", "PASS", ".", "GT", "0/1"]
        gender = "M"
        test_var = SNV(*line[:6])
        
        variant = self.vcf_loader.construct_variant(line, gender)
        
        self.assertEqual(variant.get_key(), test_var.get_key())
        self.assertFalse(hasattr(variant, "format"))
        
        # check that construct variant works for CNVs
        line = ["1", "100", ".", "T", "<DEL>", "1000", "PASS", "END=200", "GT", "0/1"]
        gender = "M"
        test_var = CNV(*line[:6])
        test_var.add_info(line[7], {})
        
        variant = self.vcf_loader.construct_variant(line, gender)
        
        self.assertEqual(variant.get_key(), test_var.get_key())
        self.assertTrue(hasattr(variant, "format"))
        
        # TODO: add checks for when HGNC is in the the filters
    
    def test_include_variant(self):
        """ check that include_variant() works correctly
        """
        
        child_variants = False
        gender = "M"
        # make a child var which passes the filters (by virtue of not having
        # initalised any filters), which should return True
        line = ["1", "100", ".", "T", "A", "1000", "PASS", ".", "GT", "0/1"]
        self.assertTrue(self.vcf_loader.include_variant(line, child_variants, gender))
        
        # make a child var that fails the filters, which should return False
        self.vcf_loader.filters = {"FILTER": ["list", "PASS", "."]} 
        line = ["1", "100", ".", "T", "A", "1000", "FAIL", ".", "GT", "0/1"]
        self.assertFalse(self.vcf_loader.include_variant(line, child_variants, gender))
        
        # now check for parents variants
        child_variants = True
        # check a parents var, where we have a matching child var
        self.vcf_loader.child_keys = set([("1", "100"), ("X", "200")])
        line = ["1", "100", ".", "T", "A", "1000", "FAIL", ".", "GT", "0/1"]
        self.assertTrue(self.vcf_loader.include_variant(line, child_variants, gender))
        
        # check a parents var, where we don't have a matching child var
        line = ["1", "200", ".", "T", "A", "1000", "FAIL", ".", "GT", "0/1"]
        self.assertFalse(self.vcf_loader.include_variant(line, child_variants, gender))
        
        # and check parental CNVs
        line = ["1", "100", ".", "T", "<DEL>", "1000", "PASS", "END=200", "GT", "0/1"]
        gender = "M"
        test_var = CNV(*line[:6])
        test_var.add_info(line[7], {})
        
        # in this function we look for overlap in CNVs. Set up a child CNV 
        # that the parents CNV must match.
        self.vcf_loader.cnv_matcher = MatchCNVs([test_var])
        self.assertTrue(self.vcf_loader.include_variant(line, child_variants, gender))
        
        # check that a parental CNV without any overlap to any childs CNVs, 
        # fails to pass
        line = ["1", "300", ".", "T", "<DEL>", "1000", "PASS", "END=400", "GT", "0/1"]
        gender = "M"
        self.assertFalse(self.vcf_loader.include_variant(line, child_variants, gender))
    
    def test_filter_de_novos(self):
        """ check that filter_de_novos() works correctly
        """
        
        # make a family without parents
        family = Family("fam_id")
        child_gender = "female"
        family.add_child("child_id", "child_vcf_path", "2", child_gender)
        self.vcf_loader.family = family
        
        # set up an autosomal variant
        line = ["1", "100", ".", "T", "G", "1000", "PASS", ".", "GT", "0/1"]
        gender = "M"
        child_var = SNV(*line[:6])
        child_var.add_info(line[7], self.vcf_loader.tags_dict)
        child_var.add_format(line[8], line[9])
        child_var.set_gender(child_gender)
        child_var.set_genotype()
        
        # combine the variant into a list of TrioGenotypes
        child_vars = [child_var]
        mother_vars = []
        father_vars = []
        trio_variants = self.vcf_loader.combine_trio_variants(child_vars, mother_vars, father_vars)
        
        # check that vars without parents get passed through automatically
        self.assertEqual(self.vcf_loader.filter_de_novos(trio_variants), trio_variants)
        
        # now add parents to the family
        family.add_mother("mother_id", "mother_vcf_path", "1", "female")
        family.add_father("father_id", "father_vcf_path", "1", "male")
        
        # re-generate the variants list now that parents have been included
        trio_variants = self.vcf_loader.combine_trio_variants(child_vars, mother_vars, father_vars)
        
        # check that vars with parents, and that appear to be de novo are
        # filtered out
        self.assertEqual(self.vcf_loader.filter_de_novos(trio_variants), [])
        
        # check that vars with parents, but which are not de novo, are retained
        mother_vars = child_vars
        trio_variants = self.vcf_loader.combine_trio_variants(child_vars, mother_vars, father_vars)
        self.assertEqual(self.vcf_loader.filter_de_novos(trio_variants), trio_variants)
        
        
        



