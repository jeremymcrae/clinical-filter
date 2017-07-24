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
import hashlib
import tempfile

from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.utils import open_vcf, get_vcf_header, exclude_header, \
    construct_variant, get_vcf_provenance
from clinicalfilter.ped import Family, Person

IS_PYTHON3 = sys.version_info.major == 3

from tests.utils import make_minimal_vcf
from tests.utils import write_temp_vcf, write_gzipped_vcf

class TestUtilsPy(unittest.TestCase):
    """ test that the utils functions work as expected
    """
    
    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.mkdtemp()
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_dir)
    
    def test_open_vcf(self):
        """ test obtaining a file handle for the VCF
        """
        
        vcf = make_minimal_vcf()
        path = os.path.join(self.temp_dir, "temp.vcf")
        write_temp_vcf(path, vcf)
        
        # check that plain VCF files can be loaded
        handle = open_vcf(path)
        self.assertEqual(type(handle), io.TextIOWrapper)
        handle.close()
        
        # check that gzipped vcf files are handled correctly
        path = os.path.join(self.temp_dir, "temp.vcf.gz")
        write_gzipped_vcf(path, vcf)
        
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
        write_temp_vcf(path, vcf)
        with self.assertRaises(OSError):
            open_vcf(path)
    
    def test_get_vcf_header(self):
        """ test that get_vcf_header() works correctly
        """
        
        vcf = make_minimal_vcf()
        path = os.path.join(self.temp_dir, "temp.vcf")
        write_temp_vcf(path, vcf)
        
        header = get_vcf_header(path)
        
        # check that the header is returned correctly
        self.assertEqual(header, vcf[:4])
    
    def test_exclude_header(self):
        """ test that exclude_header() works correctly
        """
        
        vcf = make_minimal_vcf()
        
        # make sure we drop the header, and only the header from the file
        # check this by reading the file, and making sure the first line
        # is the line we expect from the VCF
        path = os.path.join(self.temp_dir, "temp.vcf")
        write_temp_vcf(path, vcf)
        handler = open(path, "r")
        exclude_header(handler)
        self.assertEqual(handler.readline(), vcf[4])
        handler.close()
        
        # also check for gzipped VCF files.
        path = os.path.join(self.temp_dir, "temp.vcf.gz")
        write_gzipped_vcf(path, vcf)
        
        mode = 'r'
        if IS_PYTHON3:
            mode = 'rt'
        
        with gzip.open(path, mode) as handler:
            exclude_header(handler)
            self.assertEqual(handler.readline(), vcf[4])
    
    def test_get_vcf_provenance(self):
        """ test that get_vcf_provenance() works correctly
        """
        
        path = os.path.join(self.temp_dir, "temp.vcf")
        gz_path = os.path.join(self.temp_dir, "temp.vcf.gz")
        date_path = os.path.join(self.temp_dir, "temp.process.2014-02-20.vcf")
        
        family = Family('famid')
        family.add_child('child_id', 'mother', 'father', 'f', '2', path)
        family.add_mother('mom_id', '0', '0', 'female', '1', gz_path)
        family.add_father('dad_id', '0', '0', 'male', '1', date_path)
        family.set_child()
        
        vcf = make_minimal_vcf()
        vcf_string = "".join(vcf)
        if IS_PYTHON3:
            vcf_string = vcf_string.encode("utf-8")
        ungzipped_hash = hashlib.sha1(vcf_string).hexdigest()
        header = vcf[:4]
        
        write_temp_vcf(path, vcf)
        
        # check that the file defs return correctly
        (checksum, basename, date) = get_vcf_provenance(family.child)
        
        self.assertEqual(checksum, ungzipped_hash)
        self.assertEqual(basename, "temp.vcf")
        self.assertEqual(date, "2014-01-01")
        
        # now write a gzip file, and check that we get the correct hash
        write_gzipped_vcf(gz_path, vcf)
        handle = open(gz_path, "rb")
        gzipped_hash = hashlib.sha1(handle.read()).hexdigest()
        handle.close()
        
        (checksum, basename, date) = get_vcf_provenance(family.mother)
        self.assertEqual(checksum, gzipped_hash)
        
        # check that when a fileDate isn't available in the VCF, we can pick
        # the date from the path
        vcf.pop(1)
        write_temp_vcf(date_path, vcf)
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
        test_var = SNV(*line, gender=gender)
        
        variant = construct_variant(line, gender)
        
        self.assertEqual(variant.get_key(), test_var.get_key())
        self.assertEqual(variant.format, {'GT': '0/1'})
        
        # check that construct variant works for CNVs
        line = ["1", "100", ".", "T", "<DEL>", "1000", "PASS", "END=200", "GT", "0/1"]
        gender = "M"
        test_var = CNV(*line, gender=gender)
        
        variant = construct_variant(line, gender)
        
        self.assertEqual(variant.get_key(), test_var.get_key())
        self.assertEqual(variant.format, {'GT': '0/1'})
