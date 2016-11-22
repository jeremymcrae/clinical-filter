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
import copy
import shutil
import tempfile
import random
import hashlib
import subprocess

from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.trio_genotypes import TrioGenotypes
from clinicalfilter.load_vcfs import LoadVCFs
from clinicalfilter.utils import open_vcf, get_vcf_header, exclude_header, \
    construct_variant, get_vcf_provenance
from clinicalfilter.ped import Family, Person

IS_PYTHON3 = sys.version_info.major == 3

from tests.utils import make_vcf_line, make_vcf_header, make_minimal_vcf
from tests.utils import create_snv, create_cnv

class TestLoadVCFsPy(unittest.TestCase):
    """ test that the LoadVCFs methods work as expected
    """
    
    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.mkdtemp()
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_dir)
    
    def setUp(self):
        """ define a default LoadVCFs object
        """
        
        total_trios = 1
        maf_tags = ["AFR_AF", "AMR_AF", "ASN_AF", "DDD_AF", "EAS_AF", "ESP_AF",
            "EUR_AF", "MAX_AF", "SAS_AF", "UK10K_cohort_AF"]
        self.known_genes = {"ATRX": {"inheritance": {"Hemizygous": \
            {"Loss of function"}}, "start": 1, "chrom": "1", \
            "confirmed_status": {"confirmed dd gene"}, "end": 20000000}}
        
        self.vcf_loader = LoadVCFs(total_trios, maf_tags, self.known_genes, set(), None, None, )
    
    def write_temp_vcf(self, path, vcf_data):
        """ writes data to a file
        """
        
        with open(path, 'w') as handle:
            handle.writelines(vcf_data)
    
    def write_gzipped_vcf(self, path, lines):
        ''' write, compress, and index lines for a VCF
        '''
    
        with tempfile.NamedTemporaryFile(dir=self.temp_dir) as handle:
            for x in lines:
                handle.write(x.encode('utf8'))
            handle.flush()
    
            # assume bgzip and tabix binaries are available, this should be
            # handled by travis-ci setup.
            with open(path, 'w') as output:
                subprocess.call(['bgzip', '-c', handle.name], stdout=output)
            subprocess.call(['tabix', '-f', '-p', 'vcf', path])
    
    def test_open_vcf(self):
        """ test obtaining a file handle for the VCF
        """
        
        vcf = make_minimal_vcf()
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
        
        vcf = make_minimal_vcf()
        path = os.path.join(self.temp_dir, "temp.vcf")
        self.write_temp_vcf(path, vcf)
        
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
        
        vcf = make_vcf_header()
        vcf.append(make_vcf_line(pos=1, extra='HGNC=TEST;MAX_AF=0.0001'))
        vcf.append(make_vcf_line(pos=2, extra='HGNC=ATRX;MAX_AF=0.0001'))
        
        path = os.path.join(self.temp_dir, "temp.vcf")
        self.write_temp_vcf(path, vcf)
        
        person = Person('fam_id', 'sample', 'dad', 'mom', 'F', '2', path)
        
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
    
    def test_open_individual_with_mnvs(self):
        ''' test that open_individual works with MNVs
        '''
        
        vcf = make_vcf_header()
        vcf.append(make_vcf_line(pos=1, cq='splice_region_variant',
            extra='HGNC=ATRX;MAX_AF=0.0001'))
        vcf.append(make_vcf_line(pos=2, cq='missense_variant',
            extra='HGNC=ATRX;MAX_AF=0.0001'))
        
        path = os.path.join(self.temp_dir, "temp.vcf.gz")
        self.write_gzipped_vcf(path, vcf)
        
        person = Person('fam_id', 'sample', 'dad', 'mom', 'F', '2', path)
        
        args = {'chrom': "1", 'position': 1, 'id': ".", 'ref': "G", 'alts': "T",
            'filter': "PASS", 'info': "CQ=splice_region_variant;HGNC=ATRX;MAX_AF=0.0001",
            'format': "DP:GT", 'sample': "50:0/1", 'gender': "female",
            'mnv_code': 'modified_protein_altering_mnv'}
        var1 = SNV(**args)
        
        args['position'] = 2
        args['mnv_code'] = None
        args['info'] = "CQ=missense_variant;HGNC=ATRX;MAX_AF=0.0001"
        var2 = SNV(**args)
        
        # by default only one variant passes
        self.assertEqual(self.vcf_loader.open_individual(person), [var2])
        
        # if we include MNVs, then the passing variants swap
        self.assertEqual(self.vcf_loader.open_individual(person,
            mnvs={('1', 1): 'modified_protein_altering_mnv',
            ('1', 2): 'modified_synonymous_mnv'}), [var1])
    
    def test_load_trio(self):
        ''' test that load_trio() works correctly
        '''
        
        def make_vcf(person):
            # make a VCF, where one line would pass the default filtering
            vcf = make_vcf_header()
            vcf.append(make_vcf_line(pos=1, extra='HGNC=TEST;MAX_AF=0.0001'))
            vcf.append(make_vcf_line(pos=2, extra='HGNC=ATRX;MAX_AF=0.0001'))
            
            path = os.path.join(self.temp_dir, "{}.vcf.gz".format(person))
            self.write_gzipped_vcf(path, vcf)
            return path
        
        child_path = make_vcf('child')
        mother_path = make_vcf('mother')
        father_path = make_vcf('father')
        
        family = Family('fam_id')
        family.add_child('sample', 'mother_id', 'father_id', 'female', '2', child_path)
        family.add_mother('mother_id', '0', '0', 'female', '1', mother_path)
        family.add_father('father_id', '0', '0', 'male', '1', father_path)
        family.set_child()
        
        # define the parameters and values for the SNV class
        args = {'chrom': "1", 'position': 2, 'id': ".", 'ref': "G", 'alts': "T",
            'filter': "PASS", 'info': "CQ=missense_variant;HGNC=ATRX;MAX_AF=0.0001",
            'format': "DP:GT", 'sample': "50:0/1", 'gender': "female",
            'mnv_code': None}
        dad_args = copy.deepcopy(args)
        dad_args['gender'] = 'male'
        
        self.assertEqual(self.vcf_loader.load_trio(family),
            [TrioGenotypes(chrom="1", pos=2, child=SNV(**args),
                mother=SNV(**args), father=SNV(**dad_args)) ])
    
    def test_get_parental_var_snv(self):
        ''' check that get_parental_var() works correctly for SNVs
        '''
        
        sex = 'F'
        var = create_snv(sex, '0/1')
        mom = Person('fam_id', 'mom', '0', '0', 'F', '1', '/PATH')
        parental = []
        
        # try to get a matching variant for a mother. This will create a default
        # variant for a missing parental genotype
        self.assertEqual(self.vcf_loader.get_parental_var(var, parental, mom),
            SNV(chrom="1", position=150, id=".", ref="A", alts="G",
                filter="PASS", info=var.get_info_as_string(), format="GT", sample="0/0",
                gender="female", mnv_code=None))
        
        # now see if we can pick up a  variant where it does exist
        mother_var = create_snv(sex, '0/0')
        self.assertEqual(self.vcf_loader.get_parental_var(var, [mother_var],
            mom), mother_var)
    
    def test_get_parental_var_cnv(self):
        ''' check that get_parental_var() works correctly for CNVs
        '''
        
        sex = 'F'
        var = create_cnv(sex, 'deNovo')
        mom = Person('fam_id', 'mom', '0', '0', 'F', '1', '/PATH')
        parental_vars = []
        
        self.assertEqual(self.vcf_loader.get_parental_var(var, parental_vars,
            mom), CNV(chrom="1", position=150, id=".", ref="A",
                alts="<REF>", filter="PASS", info=var.get_info_as_string(), format=None,
                sample=None, gender="female", mnv_code=None))
        
        # check that even if a CNV exist in the parent at a matching site, we
        # still create a new CNV objectr for the parent
        mother_var = create_cnv(sex, 'uncertain')
        self.assertEqual(self.vcf_loader.get_parental_var(var, [mother_var],
            mom), CNV(chrom="1", position=150, id=".", ref="A",
                alts="<REF>", filter="PASS", info=var.get_info_as_string(), format=None,
                sample=None, gender="female", mnv_code=None))
    
    def test_get_parental_var_cnv_maternally_inherited(self):
        '''
        '''
        
        sex = 'F'
        mom = Person('fam_id', 'mom', '0', '0', 'F', '1', '/PATH')
        
        # check that even if a CNV exist in the parent at a matching site, we
        # still create a new CNV object for the parent
        var = create_cnv(sex, 'maternal')
        self.assertEqual(self.vcf_loader.get_parental_var(var, [], mom),
            CNV(chrom="1", position=150, id=".", ref="A",
                alts="<DUP>", filter="PASS", info=var.get_info_as_string(), format=None,
                sample=None, gender="female", mnv_code=None))
    
    def test_filter_de_novos(self):
        """ check that filter_de_novos() works correctly
        """
        
        # make a family without parents
        family = Family("fam_id")
        child_gender = "female"
        family.add_child('child_id', 'mother_id', 'father_id', child_gender, '2', 'child_path')
        self.vcf_loader.family = family
        
        # set up an autosomal variant
        gender = "M"
        args = ["1", "100", ".", "T", "G", "PASS", ".", "GT", "0/1", gender]
        child_var = SNV(*args)
        
        # combine the variant into a list of TrioGenotypes
        child_vars = [child_var]
        mother_vars = []
        father_vars = []
        trio_variants = self.vcf_loader.combine_trio_variants(family, child_vars, mother_vars, father_vars)
        
        # check that vars without parents get passed through automatically
        self.assertEqual(self.vcf_loader.filter_de_novos(trio_variants, 0.9), trio_variants)
        
        # now add parents to the family
        family.add_mother("mother_id", '0', '0', 'female', '1', "mother_vcf_path")
        family.add_father("father_id", '0', '0', 'male', '1', "father_vcf_path")
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
        
        total_trios = 1
        known_genes = {}
        maf_tags = None
        
        # if the debug info isn't available, then the SNV object doesn't use the
        # debug filter function
        self.vcf_loader = LoadVCFs(total_trios, maf_tags, known_genes, set(), None, None)
        self.assertNotEqual(SNV.passes_filters, SNV.passes_filters_with_debug)
        
        # if the debug info is passed in, check that the debug filter function
        # got set correctly
        self.vcf_loader = LoadVCFs(total_trios, maf_tags, known_genes, set(), "1", "10000")
        self.assertEqual(SNV.passes_filters, SNV.passes_filters_with_debug)
