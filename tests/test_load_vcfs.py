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
import gzip
import os
import sys
import io
import shutil
import tempfile
import random
import hashlib

from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.trio_genotypes import TrioGenotypes
from clinicalfilter.match_cnvs import MatchCNVs
from clinicalfilter.load_vcfs import LoadVCFs
from clinicalfilter.utils import open_vcf, get_vcf_header, exclude_header, \
    construct_variant, get_vcf_provenance
from clinicalfilter.ped import Family, Person

IS_PYTHON3 = sys.version_info.major == 3

class TestLoadVCFsPy(unittest.TestCase):
    """ test that the LoadVCFs methods work as expected
    """
    
    def setUp(self):
        """ define a default LoadVCFs object
        """
        
        total_trios = 1
        self.known_genes = {"ATRX": {"inheritance": {"Hemizygous": \
            {"Loss of function"}}, "start": 1, "chrom": "1", \
            "confirmed_status": {"Confirmed DD Gene"}, "end": 20000000}}
        
        self.vcf_loader = LoadVCFs(total_trios, self.known_genes, set(), None, None, )
        
        # make a temp directory for the cache file
        self.temp_dir = tempfile.mkdtemp()
    
    def tearDown(self):
        """ remove the temp directory once a test completes
        """
        
        shutil.rmtree(self.temp_dir)
    
    def make_vcf_header(self):
    
        # generate a test VCF
        lines = ['##fileformat=VCFv4.1\n',
            "##fileDate=2014-01-01\n",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
            '#CHROM\tPOS\t ID\tREF\t ALT\t QUAL\tFILTER\tINFO\tFORMAT\tsample\n']
        
        return lines
    
    def make_vcf_line(self, chrom=1, pos=1, ref='G', alts='T',
            cq='missense_variant', extra=None):
        ''' generate a VCF line suitable for the unit tests
        
        Args:
            chrom: chromosome as string
            pos: nucleotide position of the variant
            ref: reference allele
            alts: comma-separated alternate alleles
            cq: vep consequence string. Can be '|' separated (for multiple
                genes) and/or ',' separated (for multiple alt alleles).
        
        Returns:
            string for VCF line
        '''
        
        info = 'CQ={}'.format(cq)
        if extra is not None:
             info += ';' + extra
        
        return '{}\t{}\t.\t{}\t{}\t1000\tPASS\t{}\tGT:DP\t0/1:50\n'.format(chrom,
            pos, ref, alts, info)
    
    def make_minimal_vcf(self):
        """ construct the bare minimum of lines for a VCF file
        """
        
        variants = []
        variants.append(self.make_vcf_line(pos=100))
        variants.append(self.make_vcf_line(pos=200))
        
        return self.make_vcf_header() + variants
    
    def write_temp_vcf(self, path, vcf_data):
        """ writes data to a file
        """
        
        with open(path, 'w') as handle:
            handle.write("".join(vcf_data))
    
    def write_gzipped_vcf(self, path, vcf_data):
        """ writes data to a gzip file
        """
        
        mode = 'wb'
        if IS_PYTHON3:
            mode = 'wt'
        
        with gzip.open(path, mode) as handle:
            handle.write("".join(vcf_data))
    
    def test_open_vcf(self):
        """ test obtaining a file handle for the VCF
        """
        
        vcf = self.make_minimal_vcf()
        path = os.path.join(self.temp_dir, "temp.vcf")
        self.write_temp_vcf(path, vcf)
        
        # check that plain VCF files can be loaded
        handle = open_vcf(path)
        self.assertEqual(type(handle), io.TextIOWrapper)
        handle.close()
        
        # check that gzipped vcf files are handled correctly
        path = os.path.join(self.temp_dir, "temp.vcf.gz")
        self.write_gzipped_vcf(path, vcf)
        
        handle = open_vcf(path)
        if IS_PYTHON3:
            self.assertEqual(type(handle), io.TextIOWrapper)
        else:
            self.assertEqual(type(handle), gzip.GzipFile)
        handle.close()
        
        # make sure files that don't exists raise an error
        path = os.path.join(self.temp_dir, "zzz.txt")
        with self.assertRaises(OSError):
            open_vcf(path)
        
        # check that files with unknown extensions raise errors
        path = os.path.join(self.temp_dir, "temp.zzz")
        self.write_temp_vcf(path, vcf)
        with self.assertRaises(OSError):
            open_vcf(path)
    
    def test_get_vcf_header(self):
        """ test that get_vcf_header() works correctly
        """
        
        vcf = self.make_minimal_vcf()
        path = os.path.join(self.temp_dir, "temp.vcf")
        self.write_temp_vcf(path, vcf)
        
        header = get_vcf_header(path)
        
        # check that the header is returned correctly
        self.assertEqual(header, vcf[:4])
    
    def test_exclude_header(self):
        """ test that exclude_header() works correctly
        """
        
        vcf = self.make_minimal_vcf()
        
        # make sure we drop the header, and only the header from the file
        # check this by reading the file, and making sure the first line
        # is the line we expect from the VCF
        path = os.path.join(self.temp_dir, "temp.vcf")
        self.write_temp_vcf(path, vcf)
        handler = open(path, "r")
        exclude_header(handler)
        self.assertEqual(handler.readline(), vcf[4])
        handler.close()
        
        # also check for gzipped VCF files.
        path = os.path.join(self.temp_dir, "temp.vcf.gz")
        self.write_gzipped_vcf(path, vcf)
        
        mode = 'r'
        if IS_PYTHON3:
            mode = 'rt'
        
        with gzip.open(path, mode) as handler:
            exclude_header(handler)
            self.assertEqual(handler.readline(), vcf[4])
    
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
        
        path = os.path.join(self.temp_dir, "temp.vcf")
        gz_path = os.path.join(self.temp_dir, "temp.vcf.gz")
        date_path = os.path.join(self.temp_dir, "temp.process.2014-02-20.vcf")
        
        family = Family('famid')
        family.add_child('child_id', path, '2', 'F')
        family.add_mother('mom_id', gz_path, '1', 'F')
        family.add_father('mom_id', date_path, '1', 'M')
        family.set_child()
        
        vcf = self.make_minimal_vcf()
        vcf_string = "".join(vcf)
        if IS_PYTHON3:
            vcf_string = vcf_string.encode("utf-8")
        ungzipped_hash = hashlib.sha1(vcf_string).hexdigest()
        header = vcf[:4]
        
        self.write_temp_vcf(path, vcf)
        
        # check that the file defs return correctly
        (checksum, basename, date) = get_vcf_provenance(family.child)
        
        self.assertEqual(checksum, ungzipped_hash)
        self.assertEqual(basename, "temp.vcf")
        self.assertEqual(date, "2014-01-01")
        
        # now write a gzip file, and check that we get the correct hash
        self.write_gzipped_vcf(gz_path, vcf)
        handle = open(gz_path, "rb")
        gzipped_hash = hashlib.sha1(handle.read()).hexdigest()
        handle.close()
        
        (checksum, basename, date) = get_vcf_provenance(family.mother)
        self.assertEqual(checksum, gzipped_hash)
        
        # check that when a fileDate isn't available in the VCF, we can pick
        # the date from the path
        vcf.pop(1)
        self.write_temp_vcf(date_path, vcf)
        (checksum, basename, date) = get_vcf_provenance(family.father)
        self.assertEqual(date, "2014-02-20")
        
        # and check we get null values if the family member is not present
        family.father = None
        provenance = get_vcf_provenance(family.father)
        self.assertEqual(provenance, ('NA', 'NA', 'NA'))
    
    def test_construct_variant(self):
        """ test that construct_variant() works correctly
        """
        
        # check that construct variant works for SNVs
        line = ["1", "100", ".", "T", "G", "1000", "PASS", ".", "GT", "0/1"]
        gender = "M"
        test_var = SNV(*line[:6])
        
        variant = construct_variant(line, gender, self.known_genes)
        
        self.assertEqual(variant.get_key(), test_var.get_key())
        # initally constructing a SNV shouldn't affect the format variable
        self.assertEqual(variant.format, None)
        
        # check that construct variant works for CNVs
        line = ["1", "100", ".", "T", "<DEL>", "1000", "PASS", "END=200", "GT", "0/1"]
        gender = "M"
        test_var = CNV(*line[:6])
        test_var.add_info(line[7])
        
        variant = construct_variant(line, gender, self.known_genes)
        
        self.assertEqual(variant.get_key(), test_var.get_key())
        self.assertNotEqual(variant.format, None)
        
        # TODO: add checks for when HGNC is in the the filters
    
    def test_include_variant(self):
        """ check that include_variant() works correctly
        """
        
        mnvs = {}
        child_variants = False
        gender = "M"
        # make a child var which passes the filters
        line = ["1", "100", ".", "T", "A", "1000", "PASS", "CQ=missense_variant;HGNC=ATRX", "GT", "0/1"]
        self.assertTrue(self.vcf_loader.include_variant(line, child_variants, gender, mnvs))
        
        # make a child var that fails the filters, which should return False
        line = ["1", "100", ".", "T", "A", "1000", "FAIL", "CQ=missense_variant;HGNC=ATRX", "GT", "0/1"]
        self.assertFalse(self.vcf_loader.include_variant(line, child_variants, gender, mnvs))
        
        # now check for parents variants
        child_variants = True
        # check a parents var, where we have a matching child var
        self.vcf_loader.child_keys = set([("1", 100), ("X", 200)])
        line = ["1", "100", ".", "T", "A", "1000", "FAIL", "CQ=missense_variant;HGNC=ATRX", "GT", "0/1"]
        self.assertTrue(self.vcf_loader.include_variant(line, child_variants, gender, mnvs))
        
        # check a parents var, where we don't have a matching child var
        line = ["1", "200", ".", "T", "A", "1000", "FAIL", "CQ=missense_variant;HGNC=ATRX", "GT", "0/1"]
        self.assertFalse(self.vcf_loader.include_variant(line, child_variants, gender, mnvs))
        
        # and check parental CNVs
        line = ["1", "100", ".", "T", "<DEL>", "1000", "PASS", "END=200", "GT", "0/1"]
        gender = "M"
        test_var = CNV(*line[:6])
        test_var.add_info(line[7])
        
        # in this function we look for overlap in CNVs. Set up a child CNV
        # that the parents CNV must match.
        self.vcf_loader.cnv_matcher = MatchCNVs([test_var])
        self.assertTrue(self.vcf_loader.include_variant(line, child_variants, gender, mnvs))
        
        # check that a parental CNV without any overlap to any childs CNVs,
        # fails to pass
        line = ["1", "300", ".", "T", "<DEL>", "1000", "PASS", "END=400", "GT", "0/1"]
        gender = "M"
        self.assertFalse(self.vcf_loader.include_variant(line, child_variants, gender, mnvs))
    
    def test_open_individual(self):
        ''' test that open_individual() works correctly
        '''
        
        # missing individual returns empty list
        self.assertEqual(self.vcf_loader.open_individual(None), [])
        
        vcf = self.make_vcf_header()
        vcf.append(self.make_vcf_line(pos=1, extra='HGNC=TEST;MAX_AF=0.0001'))
        vcf.append(self.make_vcf_line(pos=2, extra='HGNC=ATRX;MAX_AF=0.0001'))
        
        path = os.path.join(self.temp_dir, "temp.vcf")
        self.write_temp_vcf(path, vcf)
        
        person = Person('sample', path, '2', 'F')
        
        var1 = SNV(chrom="1", position=1, id=".", ref="G", alts="T",
            filter="PASS", info="CQ=missense_variant;HGNC=TEST;MAX_AF=0.0001",
            format="DP:GT", sample="50:0/1", gender="female", mnv_code=None)
        var2 = SNV(chrom="1", position=2, id=".", ref="G", alts="T",
            filter="PASS", info="CQ=missense_variant;HGNC=ATRX;MAX_AF=0.0001",
            format="DP:GT", sample="50:0/1", gender="female", mnv_code=None)
        
        self.assertEqual(self.vcf_loader.open_individual(person), [var2])
        
        # define a set of variants to automatically pass, and check that these
        # variants pass.
        self.vcf_loader.child_keys = set([('1', 1), ('1', 2)])
        self.assertEqual(self.vcf_loader.open_individual(person,
            child_variants=True), [var1, var2])
    
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
        child_var.add_info(line[7])
        child_var.add_format(line[8], line[9])
        child_var.set_gender(child_gender)
        child_var.set_genotype()
        
        # combine the variant into a list of TrioGenotypes
        child_vars = [child_var]
        mother_vars = []
        father_vars = []
        trio_variants = self.vcf_loader.combine_trio_variants(family, child_vars, mother_vars, father_vars)
        
        # check that vars without parents get passed through automatically
        self.assertEqual(self.vcf_loader.filter_de_novos(trio_variants, 0.9), trio_variants)
        
        # now add parents to the family
        family.add_mother("mother_id", "mother_vcf_path", "1", "female")
        family.add_father("father_id", "father_vcf_path", "1", "male")
        self.vcf_loader.family = family
        
        # re-generate the variants list now that parents have been included
        trio_variants = self.vcf_loader.combine_trio_variants(family, child_vars, mother_vars, father_vars)
        
        # check that vars with parents, and that appear to be de novo are
        # filtered out
        self.assertEqual(self.vcf_loader.filter_de_novos(trio_variants, 0.9), [])
        
        # check that vars with parents, but which are not de novo, are retained
        mother_vars = child_vars
        trio_variants = self.vcf_loader.combine_trio_variants(family, child_vars, mother_vars, father_vars)
        
        self.assertEqual(self.vcf_loader.filter_de_novos(trio_variants, 0.9), trio_variants)
    
    def test_debug_option(self):
        """ test whether we can set up the class with the debug option
        """
        
        counter = 0
        total_trios = 1
        known_genes = {}
        
        self.vcf_loader = LoadVCFs(total_trios, known_genes, set(), "1", "10000")
        
        # check that the debug filter function got set correctly
        self.assertEqual(SNV.passes_filters, SNV.passes_filters_with_debug)
