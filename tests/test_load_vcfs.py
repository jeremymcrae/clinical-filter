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

from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.variant.info import Info
from clinicalfilter.trio_genotypes import TrioGenotypes
from clinicalfilter.load_vcfs import load_variants, include_variant, \
    open_individual, load_trio, combine_trio_variants, get_parental_var, \
    filter_de_novos
from clinicalfilter.ped import Family, Person

IS_PYTHON3 = sys.version_info.major == 3

from tests.utils import make_vcf_line, make_vcf_header, make_minimal_vcf
from tests.utils import create_snv, create_cnv, write_temp_vcf, write_gzipped_vcf

class TestLoadVCFsPy(unittest.TestCase):
    """ test that the vcf loading functions work as expected
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
        
        self.known_genes = {"ATRX": {"inheritance": {"Hemizygous": \
            {"Loss of function"}}, "start": 1, "chrom": "1", \
            "confirmed_status": {"confirmed dd gene"}, "end": 20000000}}
        
        # assign the known genes to the variant classes, so later tests can rely
        # on the known genes
        SNV.known_genes = self.known_genes
        CNV.known_genes = self.known_genes
    
    def tearDown(self):
        # assign the known genes to the variant classes, so later tests can rely
        # on the known genes
        SNV.known_genes = self.known_genes
        CNV.known_genes = self.known_genes
        
        Info.populations = []
        Info.last_base = set()
    
    def test_load_variants(self):
        ''' test that load_variants() works correctly. Mainly checks variables are set correctly
        '''
        
        vcf = make_minimal_vcf()
        path = os.path.join(self.temp_dir, "temp.vcf.gz")
        write_gzipped_vcf(path, vcf)
        
        fam = Family('fam', children=[Person('fam', 'child', '0', '0', 'f', '2', path)])
        variants = load_variants(fam, 0.9, ['AFR_AF'], self.known_genes, set())
        
        self.assertEqual(SNV.known_genes, self.known_genes)
        self.assertEqual(CNV.known_genes, self.known_genes)
        self.assertEqual(Info.populations, ['AFR_AF'])
        self.assertEqual(Info.last_base, set())
        
        # and check that the
        variants = load_variants(fam, 0.9, [], None, set([('1', 100)]))
        self.assertIsNone(SNV.known_genes, self.known_genes)
        self.assertIsNone(CNV.known_genes, self.known_genes)
        self.assertEqual(Info.populations, [])
        self.assertEqual(Info.last_base, set([('1', 100)]))
    
    def test_include_variant(self):
        """ check that include_variant() works correctly
        """
        
        mnvs = {}
        child_keys = None
        gender = "M"
        # make a child var which passes the filters
        line = ["1", "100", ".", "T", "A", "1000", "PASS", "CQ=missense_variant;HGNC=ATRX", "GT", "0/1"]
        self.assertTrue(include_variant(line, child_keys, gender, mnvs))
        
        # make a child var that fails the filters, which should return False
        line = ["1", "100", ".", "T", "A", "1000", "FAIL", "CQ=missense_variant;HGNC=ATRX", "GT", "0/1"]
        self.assertFalse(include_variant(line, child_keys, gender, mnvs))
        
        # now check for parents variants
        # check a parents var, where we have a matching child var
        child_keys = set([("1", 100), ("X", 200)])
        line = ["1", "100", ".", "T", "A", "1000", "FAIL", "CQ=missense_variant;HGNC=ATRX", "GT", "0/1"]
        self.assertTrue(include_variant(line, child_keys, gender, mnvs))
        
        # check a parents var, where we don't have a matching child var
        line = ["1", "200", ".", "T", "A", "1000", "FAIL", "CQ=missense_variant;HGNC=ATRX", "GT", "0/1"]
        self.assertFalse(include_variant(line, child_keys, gender, mnvs))
        
        # and check parental CNVs
        line = ["1", "100", ".", "T", "<DEL>", "1000", "PASS", "END=200", "GT", "0/1"]
        gender = "M"
        test_var = CNV(*line)
        
        # in this function we look for overlap in CNVs. Set up a child CNV
        # that the parents CNV must match.
        self.assertTrue(include_variant(line, child_keys, gender, mnvs))
        
        # check that a parental CNV without any overlap to any childs CNVs,
        # fails to pass
        line = ["1", "300", ".", "T", "<DEL>", "1000", "PASS", "END=400", "GT", "0/1"]
        gender = "M"
        self.assertFalse(include_variant(line, child_keys, gender, mnvs))
    
    def test_open_individual(self):
        ''' test that open_individual() works correctly
        '''
        
        # missing individual returns empty list
        self.assertEqual(open_individual(None), [])
        
        vcf = make_vcf_header()
        vcf.append(make_vcf_line(pos=1, extra='HGNC=TEST;MAX_AF=0.0001'))
        vcf.append(make_vcf_line(pos=2, extra='HGNC=ATRX;MAX_AF=0.0001'))
        
        path = os.path.join(self.temp_dir, "temp.vcf")
        write_temp_vcf(path, vcf)
        
        person = Person('fam_id', 'sample', 'dad', 'mom', 'F', '2', path)
        
        var1 = SNV(chrom="1", position=1, id=".", ref="G", alts="T",
            qual='1000', filter="PASS", info="CQ=missense_variant;HGNC=TEST;MAX_AF=0.0001",
            format="DP:GT", sample="50:0/1", gender="female", mnv_code=None)
        var2 = SNV(chrom="1", position=2, id=".", ref="G", alts="T",
            qual='1000', filter="PASS", info="CQ=missense_variant;HGNC=ATRX;MAX_AF=0.0001",
            format="DP:GT", sample="50:0/1", gender="female", mnv_code=None)
        
        self.assertEqual(open_individual(person), [var2])
        
        # define a set of variants to automatically pass, and check that these
        # variants pass.
        child_keys = set([('1', 1), ('1', 2)])
        self.assertEqual(open_individual(person,
            child_variants=child_keys), [var1, var2])
    
    def test_open_individual_with_mnvs(self):
        ''' test that open_individual works with MNVs
        '''
        
        vcf = make_vcf_header()
        vcf.append(make_vcf_line(pos=1, cq='splice_region_variant',
            extra='HGNC=ATRX;MAX_AF=0.0001'))
        vcf.append(make_vcf_line(pos=2, cq='missense_variant',
            extra='HGNC=ATRX;MAX_AF=0.0001'))
        
        path = os.path.join(self.temp_dir, "temp.vcf.gz")
        write_gzipped_vcf(path, vcf)
        
        person = Person('fam_id', 'sample', 'dad', 'mom', 'F', '2', path)
        
        args = {'chrom': "1", 'position': 1, 'id': ".", 'ref': "G", 'alts': "T",
            'filter': "PASS", 'info': "CQ=splice_region_variant;HGNC=ATRX;MAX_AF=0.0001",
            'format': "DP:GT", 'sample': "50:0/1", 'gender': "female",
            'mnv_code': 'modified_protein_altering_mnv', 'qual': '1000'}
        var1 = SNV(**args)
        
        args['position'] = 2
        args['mnv_code'] = None
        args['info'] = "CQ=missense_variant;HGNC=ATRX;MAX_AF=0.0001"
        var2 = SNV(**args)
        
        # by default only one variant passes
        self.assertEqual(open_individual(person), [var2])
        
        # if we include MNVs, then the passing variants swap
        self.assertEqual(open_individual(person,
            mnvs={('1', 1): 'modified_protein_altering_mnv',
            ('1', 2): 'modified_synonymous_mnv'}), [var1])
    
    def test_open_individual_male_het_chrx(self):
        """ test that open_individual() passes over hets in males on chrX
        """
        
        # the sub-functions are all tested elsewhere, this test merely checks
        # that valid variants are added to the variants list, and invalid
        # variants are passed over without being added to the variants list
        
        vcf = make_vcf_header()
        vcf.append(make_vcf_line(chrom='X', pos=1, genotype='0/1',
            extra='HGNC=TEST;MAX_AF=0.0001'))
        
        path = os.path.join(self.temp_dir, "temp.vcf")
        write_temp_vcf(path, vcf)
        
        person = Person('fam_id', 'sample', 'dad', 'mom', 'M', '2', path)
        
        self.assertEqual(open_individual(person), [])
    
    def test_load_trio(self):
        ''' test that load_trio() works correctly
        '''
        
        def make_vcf(person):
            # make a VCF, where one line would pass the default filtering
            vcf = make_vcf_header()
            vcf.append(make_vcf_line(pos=1, extra='HGNC=TEST;MAX_AF=0.0001'))
            vcf.append(make_vcf_line(pos=2, extra='HGNC=ATRX;MAX_AF=0.0001'))
            
            path = os.path.join(self.temp_dir, "{}.vcf.gz".format(person))
            write_gzipped_vcf(path, vcf)
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
            'mnv_code': None, 'qual': '1000'}
        dad_args = copy.deepcopy(args)
        dad_args['gender'] = 'male'
        
        self.assertEqual(load_trio(family),
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
        self.assertEqual(get_parental_var(var, parental, mom),
            SNV(chrom="1", position=150, id=".", ref="A", alts="G",
                qual='1000', filter="PASS", info=str(var.info), format="GT",
                sample="0/0", gender="female", mnv_code=None))
        
        # now see if we can pick up a  variant where it does exist
        mother_var = create_snv(sex, '0/0')
        self.assertEqual(get_parental_var(var, [mother_var],
            mom), mother_var)
    
    def test_get_parental_var_cnv(self):
        ''' check that get_parental_var() works correctly for CNVs
        '''
        
        sex = 'F'
        var = create_cnv(sex, 'deNovo')
        mom = Person('fam_id', 'mom', '0', '0', 'F', '1', '/PATH')
        parental_vars = []
        
        self.assertEqual(get_parental_var(var, parental_vars,
            mom), CNV(chrom="1", position=150, id=".", ref="A",
                alts="<REF>", qual='1000', filter="PASS", info=str(var.info),
                format='INHERITANCE', sample='uncertain', gender="female",
                mnv_code=None))
        
        # check that even if a CNV exist in the parent at a matching site, we
        # still create a new CNV objectr for the parent
        mother_var = create_cnv(sex, 'uncertain')
        self.assertEqual(get_parental_var(var, [mother_var],
            mom), CNV(chrom="1", position=150, id=".", ref="A",
                alts="<REF>", qual='1000', filter="PASS", info=str(var.info),
                format='INHERITANCE', sample='uncertain', gender="female",
                mnv_code=None))
    
    def test_get_parental_var_cnv_maternally_inherited(self):
        ''' test that we can construct a maternally inherited CNV
        '''
        
        sex = 'F'
        mom = Person('fam_id', 'mom', '0', '0', 'F', '1', '/PATH')
        
        # check that even if a CNV exist in the parent at a matching site, we
        # still create a new CNV object for the parent
        var = create_cnv(sex, 'maternal')
        self.assertEqual(get_parental_var(var, [], mom),
            CNV(chrom="1", position=150, id=".", ref="A",
                alts="<DUP>", qual='1000',filter="PASS", info=str(var.info),
                format='INHERITANCE', sample='uncertain', gender="female",
                mnv_code=None))
    
    def test_filter_de_novos(self):
        """ check that filter_de_novos() works correctly
        """
        
        # make a family without parents
        family = Family("fam_id")
        child_gender = "female"
        family.add_child('child_id', 'mother_id', 'father_id', child_gender, '2', 'child_path')
        
        # set up an autosomal variant
        gender = "M"
        args = ["1", "100", ".", "T", "G", "1000", "PASS", ".", "GT", "0/1", gender]
        child_var = SNV(*args)
        
        # combine the variant into a list of TrioGenotypes
        child_vars = [child_var]
        mother_vars = []
        father_vars = []
        trio_variants = combine_trio_variants(family, child_vars, mother_vars, father_vars)
        
        # check that vars without parents get passed through automatically
        self.assertEqual(filter_de_novos(trio_variants, 0.9), trio_variants)
        
        # now add parents to the family
        family.add_mother("mother_id", '0', '0', 'female', '1', "mother_vcf_path")
        family.add_father("father_id", '0', '0', 'male', '1', "father_vcf_path")
        family = family
        
        # re-generate the variants list now that parents have been included
        trio_variants = combine_trio_variants(family, child_vars, mother_vars, father_vars)
        
        # check that vars with parents, and that appear to be de novo are
        # filtered out
        self.assertEqual(filter_de_novos(trio_variants, 0.9), [])
        
        # check that vars with parents, but which are not de novo, are retained
        mother_vars = child_vars
        trio_variants = combine_trio_variants(family, child_vars, mother_vars, father_vars)
        
        self.assertEqual(filter_de_novos(trio_variants, 0.9), trio_variants)
    
    def test_debug_option(self):
        """ test whether we can set up the class with the debug option
        """
        
        known = {}
        pops = None
        
        vcf = make_minimal_vcf()
        path = os.path.join(self.temp_dir, "temp.vcf.gz")
        write_gzipped_vcf(path, vcf)
        
        fam = Family('fam', children=[Person('fam', 'child', '0', '0', 'f', '2', path)])
        
        # if the debug info isn't available, then the SNV object doesn't use the
        # debug filter function
        variants = load_variants(fam, 1.0, pops, known, set())
        self.assertNotEqual(SNV.passes_filters, SNV.passes_filters_with_debug)
        
        # if the debug info is passed in, check that the debug filter function
        # got set correctly
        variants = load_variants(fam, 1.0, pops, known, set(), "1", "10000")
        self.assertEqual(SNV.passes_filters, SNV.passes_filters_with_debug)
