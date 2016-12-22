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

import os
import shutil
import tempfile
import unittest
from collections import namedtuple
import subprocess
import re

import tabix

from clinicalfilter.utils import open_vcf, exclude_header
from clinicalfilter.multinucleotide_variants import get_mnv_candidates, \
    find_nearby_variants, parse_vcf_line, get_matches, is_not_indel, is_coding, \
    screen_pairs, same_aa, translate, get_codons, check_mnv_consequence

from tests.utils import make_vcf_header, make_vcf_line

class TestMNVChecksPy(unittest.TestCase):
    """ test the multinucleotide variant (MNV) checking functions
    """
    
    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.mkdtemp()
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tempdir)
    
    def setUp(self):
        self.vcf = tempfile.NamedTemporaryFile(suffix='.vcf.gz',
            dir=self.tempdir, delete=False, mode='w')
        self.path = self.vcf.name
        
        self.Variant = namedtuple('Variant', ['chrom', 'pos', 'id', 'ref',
            'alts', 'qual', 'filter', 'info'])
        
        self.pattern = re.compile('[ACGT]')
    
    def write_vcf(self, lines):
        ''' write, compress, and index lines for a VCF
        '''
        
        with tempfile.NamedTemporaryFile(dir=self.tempdir) as handle:
            for x in lines:
                handle.write(x.encode('utf8'))
            handle.flush()
            
            # assume bgzip and tabix binaries are available, this should be
            # handled by travis-ci setup.
            subprocess.call(['bgzip', '-c', handle.name], stdout=self.vcf)
            subprocess.call(['tabix', '-f', '-p', 'vcf', self.path])
    
    def make_vcf_header(self):
    
        # generate a test VCF
        lines = ['##fileformat=VCFv4.1\n',
            '#CHROM\tPOS\t ID\tREF\t ALT\t QUAL\tFILTER\tINFO\tFORMAT\tsample\n']
        
        return lines
    
    def make_vcf_line(self, chrom=1, pos=1, ref='G', alts='T', aa_pos='1',
            codons='aGt/aTt', cq='missense_variant'):
        ''' generate a VCF line suitable for the unit tests
        
        Args:
            chrom: chromosome as string
            pos: nucleotide position of the variant
            ref: reference allele
            alts: comma-separated alternate alleles
            aa_pos: amino acid position of variant in protein sequence
            codons: codon sequence for reference and alternate alleles. The
                position of the modified base is given by the base in upper
                case. The reference and alternate codons are '/' separated.
                Multiple codons can be '|' separated (for multiple genes, and
                genes where the variant is not coding use the '.' null value).
                Codons can also be ',' separated for multiple alt alleles.
            cq: vep consequence string. Can be '|' separated (for multiple
                genes) and/or ',' separated (for multiple alt alleles).
        
        Returns:
            string for VCF line
        '''
        
        aa_string = ''
        if aa_pos is not None:
            aa_string = ';Protein_position={}'.format(aa_pos)
        
        info = 'Codons={};CQ={}{}'.format(codons, cq, aa_string)
        
        return '{}\t{}\t.\t{}\t{}\t1000\tPASS\t{}\tGT:DP\t0/1:50\n'.format(chrom,
            pos, ref, alts, info)
    
    def test_get_mnv_candidates(self):
        ''' check that get_mnv_candidates works correctly
        '''
        
        lines = make_vcf_header()
        lines.append(make_vcf_line(chrom='1', pos=1, extra='Protein_position=1;Codons=aaT/aaG'))
        lines.append(make_vcf_line(chrom='1', pos=2, extra='Protein_position=1;Codons=Aat/Cat'))
        self.write_vcf(lines)
        
        self.assertEqual(get_mnv_candidates(self.path), {
            ('1', 1): 'alternate_residue_mnv', ('1', 2): 'alternate_residue_mnv'})
    
    def test_get_mnv_candidates_catch_assertion_error(self):
        ''' check that get_mnv_candidates works correctly
        '''
        
        lines = make_vcf_header()
        lines.append(make_vcf_line(chrom='1', pos=1, extra='Protein_position=1;Codons=aaT/aaG'))
        lines.append(make_vcf_line(chrom='1', pos=2, extra='Protein_position=2;Codons=Att/Ctt'))
        self.write_vcf(lines)
        
        self.assertEqual(get_mnv_candidates(self.path), {})
    
    def test_find_nearby_variants(self):
        ''' test that find_nearby_variants() works correctly
        '''
        
        lines = make_vcf_header()
        lines.append(make_vcf_line(pos=1))
        lines.append(make_vcf_line(pos=2))
        self.write_vcf(lines)
        
        vcf = open_vcf(self.path)
        exclude_header(vcf)
        self.assertEqual(find_nearby_variants(vcf), [[('1', 1), ('1', 2)]])
    
    def test_find_nearby_variants_separated(self):
        ''' test that find_nearby_variants() doesn't include vars far apart
        '''
        
        lines = make_vcf_header()
        lines.append(make_vcf_line(pos=1))
        lines.append(make_vcf_line(pos=4))
        self.write_vcf(lines)
        
        vcf = open_vcf(self.path)
        exclude_header(vcf)
        self.assertEqual(find_nearby_variants(vcf), [])
    
    def test_find_nearby_variants_duplicate_position(self):
        ''' test that find_nearby_variants() works correctly with a duplicate var
        '''
        
        # get the default two variants
        lines = make_vcf_header()
        lines.append(make_vcf_line(pos=1))
        lines.append(make_vcf_line(pos=2))
        
        # make a third variant, but at the same position as the second
        lines.append(make_vcf_line(pos=2))
        self.write_vcf(lines)
        
        vcf = open_vcf(self.path)
        exclude_header(vcf)
        self.assertEqual(find_nearby_variants(vcf), [[('1', 1), ('1', 2)]])
    
    def test_find_nearby_variants_different_chroms(self):
        ''' test that find_nearby_variants() works correctly with successive
        variants on different chroms, but at the same position.
        '''
        
        # get the default two variants
        lines = make_vcf_header()
        lines.append(make_vcf_line(chrom='1', pos=1))
        lines.append(make_vcf_line(chrom='2', pos=1))
        
        vcf = open_vcf(self.path)
        exclude_header(vcf)
        self.assertEqual(find_nearby_variants(vcf), [])
    
    def test_find_nearby_variants_different_threshold(self):
        ''' test that find_nearby_variants() works correctly when we change the threshold distance.
        '''
        
        # get the default two variants
        lines = make_vcf_header()
        lines.append(make_vcf_line(pos=1))
        lines.append(make_vcf_line(pos=2))
        
        vcf = open_vcf(self.path)
        exclude_header(vcf)
        
        # using a lower threshold shouldn't allow any of the variants to pass
        self.assertEqual(find_nearby_variants(vcf, threshold=0), [])
    
    def test_parse_vcf_line(self):
        ''' test that parse_vcf_line() works correctly
        '''
        
        line = make_vcf_line(extra='Protein_position=1;Codons=aGt/aTt').split('\t')
        var = parse_vcf_line(line, self.Variant)
        
        parsed = self.Variant(chrom='1', pos=1, id='.', ref='G', alts=['T'],
            qual='1000', filter='PASS', info={'Protein_position': '1',
            'CQ': 'missense_variant', 'Codons': 'aGt/aTt'})
        
        self.assertEqual(var, parsed)
        
        # check that passing in an unsplit string raises an error
        with self.assertRaises(ValueError):
            line = self.make_vcf_line()
            parse_vcf_line(line, self.Variant)
    
    def test_parse_vcf_line_flag(self):
        ''' check that parse_vcf_line() works when info includes a flag
        '''
        
        line = '1\t866\t.\tCT\tCCCCTC\t200\tPASS\tAC=2;AC_Adj=500;AC_Het=500;' \
            'TEST\tGT:AD\t1/1:0,4'
        line = line.split('\t')
        var = parse_vcf_line(line, self.Variant)
        
        parsed = self.Variant(chrom='1', pos=866, id='.', ref='CT',
            alts=['CCCCTC'], qual='200', filter='PASS', info={'AC': '2',
            'AC_Adj': '500', 'AC_Het': '500', 'TEST': True})
        
        self.assertEqual(var, parsed)
    
    def test_parse_vcf_line_multi_alts(self):
        ''' check that parse_vcf_line() works when we have multiple alts
        '''
        
        line = make_vcf_line(alts='C,CT').split('\t')
        var = parse_vcf_line(line, self.Variant)
        
        parsed = self.Variant(chrom='1', pos=1, id='.', ref='G', alts=['C', 'CT'],
            qual='1000', filter='PASS', info={'CQ': 'missense_variant'})
        
        self.assertEqual(var, parsed)
    
    def test_get_matches(self):
        ''' check that get_matches works correctly
        '''
        
        # get the VCF lines
        lines = make_vcf_header()
        lines.append(make_vcf_line(pos=1))
        lines.append(make_vcf_line(pos=2))
        lines.append(make_vcf_line(pos=4))
        lines.append(make_vcf_line(pos=5))
        self.write_vcf(lines)
        
        vcf = tabix.open(self.path)
        pair = [('1', 2), ('1', 4)]
        
        # define the expected lines
        var1 = parse_vcf_line(make_vcf_line(pos=2).split('\t'), self.Variant)
        var2 = parse_vcf_line(make_vcf_line(pos=4).split('\t'), self.Variant)
        
        self.assertEqual(list(get_matches(vcf, pair)), [var1, var2])
    
    def test_get_matches_extra(self):
        ''' check that get_matches works correctly with > 2 in the 'pair'
        '''
        
        # get the VCF lines
        lines = make_vcf_header()
        lines.append(make_vcf_line(pos=1))
        lines.append(make_vcf_line(pos=2))
        lines.append(make_vcf_line(pos=4))
        lines.append(make_vcf_line(pos=5))
        self.write_vcf(lines)
        
        vcf = tabix.open(self.path)
        pair = [('1', 2), ('1', 4), ('1', 5)]
        
        # define the expected lines
        var1 = parse_vcf_line(make_vcf_line(pos=2).split('\t'), self.Variant)
        var2 = parse_vcf_line(make_vcf_line(pos=4).split('\t'), self.Variant)
        var3 = parse_vcf_line(make_vcf_line(pos=5).split('\t'), self.Variant)
        
        self.assertEqual(list(get_matches(vcf, pair)), [var1, var2, var3])
    
    def test_is_not_indel(self):
        ''' check that is_not_indel() works correctly
        '''
        
        # check a deletion indel (from ref allele)
        line = make_vcf_line(ref='AA').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertFalse(is_not_indel(var))
        
        # check a SNV, should pass
        line = make_vcf_line(ref='A').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertTrue(is_not_indel(var))
        
        # check a SNV with multiple alts
        line = make_vcf_line(ref='A', alts='T,G').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertTrue(is_not_indel(var))
        
        # check a variant with multiple alts, only one of which is for a SNV
        line = make_vcf_line(ref='A', alts='TT,G').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertTrue(is_not_indel(var))
        
        # check an indel with multiple alts, none of which are for a SNV
        line = make_vcf_line(ref='A', alts='TT,*').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertFalse(is_not_indel(var))
        
        # check a deletion indel
        line = make_vcf_line(ref='A', alts='TT').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertFalse(is_not_indel(var))
    
    def test_is_coding(self):
        ''' check that is_coding() works correctly
        '''
        
        # check for a single transcript and variant in CDS
        line = make_vcf_line(cq='missense_variant').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertTrue(is_coding(var))
        
        # check for a single transcript and synonymous variant
        line = make_vcf_line(cq='synonymous_variant').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertTrue(is_coding(var))
        
        # check for a single transcript and coding variant
        line = make_vcf_line(cq='intergenic_variant').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertFalse(is_coding(var))
        
        # check for a single transcript and coding variant
        line = make_vcf_line(cq='intergenic_variant|missense_variant').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertTrue(is_coding(var))
        
        # check for a single transcript and coding variant
        line = make_vcf_line(cq='intergenic_variant,missense_variant').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertTrue(is_coding(var))
        
        # check for a single transcript and coding variant
        line = make_vcf_line(cq='intergenic_variant,'
            'intergenic_variant|missense_variant').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertTrue(is_coding(var))
        
        # check for a single transcript and coding variant
        line = make_vcf_line(cq='intergenic_variant,'
            'intergenic_variant|intergenic_variant').split('\t')
        var = parse_vcf_line(line, self.Variant)
        self.assertFalse(is_coding(var))
    
    def test_screen_pairs(self):
        ''' test that screen_pairs() works correctly
        '''
        
        # get the VCF lines
        lines = make_vcf_header()
        lines.append(make_vcf_line(pos=1))
        lines.append(make_vcf_line(pos=2))
        lines.append(make_vcf_line(pos=4))
        lines.append(make_vcf_line(pos=5))
        lines.append(make_vcf_line(pos=7))
        lines.append(make_vcf_line(pos=8))
        self.write_vcf(lines)
        
        vcf = tabix.open(self.path)
        pairs = [[('1', 2), ('1', 4)], [('1', 7), ('1', 8)]]
        
        self.assertEqual(screen_pairs(vcf, pairs, is_not_indel), pairs)
        
        # check that the other filter function also works cleanly
        self.assertEqual(screen_pairs(vcf, pairs, is_coding), pairs)
    
    def test_screen_pairs_nonstandard_pair(self):
        ''' test that screen_pairs() works correctly
        '''
        
        # get the VCF lines
        lines = make_vcf_header()
        lines.append(make_vcf_line(pos=2))
        lines.append(make_vcf_line(pos=4))
        lines.append(make_vcf_line(pos=5))
        lines.append(make_vcf_line(pos=7))
        lines.append(make_vcf_line(pos=8))
        self.write_vcf(lines)
        
        vcf = tabix.open(self.path)
        # set up a list of 'pairs', where one 'pair' has three variants in it.
        # we exclude 'pairs' where n != 2.
        pairs = [[('1', 2), ('1', 4), ('1', 5)], [('1', 7), ('1', 8)]]
        self.assertEqual(screen_pairs(vcf, pairs, is_not_indel), [[('1', 7), ('1', 8)]])
    
    def test_same_aa(self):
        ''' check that same_aa() works correctly
        '''
        
        # get the VCF lines
        lines = make_vcf_header()
        lines.append(make_vcf_line(pos=2, extra='Protein_position=1'))
        lines.append(make_vcf_line(pos=4, extra='Protein_position=1'))
        self.write_vcf(lines)
        
        vcf = tabix.open(self.path)
        pairs = [[('1', 2), ('1', 4)]]
        
        self.assertEqual(same_aa(vcf, pairs), [[('1', 2), ('1', 4)]])
    
    def test_same_aa_different_positions(self):
        ''' check that same_aa() works correctly for different amino acids
        '''
        
        lines = make_vcf_header()
        lines.append(make_vcf_line(pos=5, extra='Protein_position=2'))
        lines.append(make_vcf_line(pos=7, extra='Protein_position=3'))
        lines.append(make_vcf_line(pos=8, extra='Protein_position=4'))
        self.write_vcf(lines)
        
        vcf = tabix.open(self.path)
        pairs = [[('1', 7), ('1', 8)]]
        
        self.assertEqual(same_aa(vcf, pairs), [])
    
    def test_same_aa_missing_protein_positions(self):
        ''' check that same_aa() works correctly when the vars aren't in the CDS
        '''
        
        # if one of the variants in the pair does not have a protein position
        # listed (i.e. residue number), that indicates the variant could be
        # affecting the splice site, so we can't use the pair.
        lines = make_vcf_header()
        lines.append(make_vcf_line(pos=5))
        lines.append(make_vcf_line(pos=7))
        lines.append(make_vcf_line(pos=8, extra='Protein_position=4'))
        self.write_vcf(lines)
        
        vcf = tabix.open(self.path)
        pairs = [[('1', 7), ('1', 8)]]
        
        self.assertEqual(same_aa(vcf, pairs), [])
    
    def test_translate(self):
        """ test that translate() works correctly
        """
        
        self.assertEqual(translate('AAG'), 'K')
        self.assertEqual(translate('TAG'), '*')
        
        # raise errors for unknown sequences and short codons
        with self.assertRaises(KeyError):
            translate('ZZZ')
            translate('AA')
    
    def test_get_codons(self):
        ''' test that get_codons() works correctly
        '''
        
        var1 = make_vcf_line(extra='Codons=aGt/aTt').split('\t')
        var2 = make_vcf_line(extra='Codons=agT/agC').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        self.assertEqual(get_codons(var1, var2, self.pattern),
            {'reference': 'aGt', 'snv1': 'aTt', 'snv2': 'agC', 'mnv': 'aTC'})
    
    def test_get_codons_null_value(self):
        ''' test that get_codons() works correctly when some genes have null values
        '''
        
        var1 = make_vcf_line(extra='Codons=aGt/aTt|.').split('\t')
        var2 = make_vcf_line(extra='Codons=agT/agC|.').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        self.assertEqual(get_codons(var1, var2, self.pattern),
            {'reference': 'aGt', 'snv1': 'aTt', 'snv2': 'agC', 'mnv': 'aTC'})
    
    def test_get_codons_duplicate_codons(self):
        ''' test that get_codons() works when variants duplicate codons
        '''
        
        var1 = make_vcf_line(extra='Codons=aGt/aTt|aGt/aTt').split('\t')
        var2 = make_vcf_line(extra='Codons=agT/agC|agT/agC').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        self.assertEqual(get_codons(var1, var2, self.pattern),
            {'reference': 'aGt', 'snv1': 'aTt', 'snv2': 'agC', 'mnv': 'aTC'})
        
    def test_get_codons_with_different_codons(self):
        ''' test that get_codons() raises an error when different transcripts
        have different codons
        '''
        
        var1 = make_vcf_line(extra='Codons=aGt/aTt|Gta/Tta').split('\t')
        var2 = make_vcf_line(extra='Codons=agT/agC|gTa/gCa').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        with self.assertRaises(AssertionError):
            get_codons(var1, var2, self.pattern)
    
    def test_get_codons_without_uppercase_base(self):
        ''' test that get_codons() raises an error when the variant position is
        not un upper case.
        '''
        
        var1 = make_vcf_line(extra='Codons=agt/att').split('\t')
        var2 = make_vcf_line(extra='Codons=agt/agc').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        with self.assertRaises(AttributeError):
            get_codons(var1, var2, self.pattern)
    
    def test_get_codons_short_codon(self):
        ''' test that get_codons() raises an error the codons are not 3bp long.
        '''
        
        var1 = make_vcf_line(extra='Codons=aG/aT').split('\t')
        var2 = make_vcf_line(extra='Codons=Ag/Gg').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        with self.assertRaises(AssertionError):
            get_codons(var1, var2, self.pattern)
    
    def test_check_minv_consequence_unmodified_synonymous(self):
        ''' test that get_mnv_consequence() works correctly
        '''
        
        var1 = make_vcf_line(extra='Codons=Cga/Aga').split('\t')
        var2 = make_vcf_line(extra='Codons=cgA/cgG').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        self.assertEqual(check_mnv_consequence(var1, var2, self.pattern),
            'unmodified_synonymous_mnv')
    
    def test_check_mnv_consequence_unmodified_altering(self):
        ''' test that get_mnv_consequence() works correctly
        '''
        
        var1 = make_vcf_line(extra='Codons=Ctt/Ttt').split('\t')
        var2 = make_vcf_line(extra='Codons=ctT/ctC').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        self.assertEqual(check_mnv_consequence(var1, var2, self.pattern),
            'unmodified_protein_altering_mnv')
        
    def test_check_mnv_consequence_modified_altering(self):
        ''' test that get_mnv_consequence() works correctly
        '''
        
        var1 = make_vcf_line(extra='Codons=Cta/Tta').split('\t')
        var2 = make_vcf_line(extra='Codons=ctA/ctT').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        self.assertEqual(check_mnv_consequence(var1, var2, self.pattern),
            'modified_protein_altering_mnv')
        
    def test_check_mnv_consequence_modified_synonymous(self):
        ''' test that get_mnv_consequence() works correctly
        
        This should only be true for Serine residues, such as TCT -> AGT.
        '''
        
        var1 = make_vcf_line(extra='Codons=tCt/tGt').split('\t')
        var2 = make_vcf_line(extra='Codons=Tct/Act').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        self.assertEqual(check_mnv_consequence(var1, var2, self.pattern),
            'modified_synonymous_mnv')
        
    def test_check_mnv_consequence_modified_stop_gained(self):
        ''' test that get_mnv_consequence() works correctly
        '''
        
        var1 = make_vcf_line(extra='Codons=Cat/Tat').split('\t')
        var2 = make_vcf_line(extra='Codons=caT/caG').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        self.assertEqual(check_mnv_consequence(var1, var2, self.pattern),
            'modified_stop_gained_mnv')
    
    def test_check_mnv_consequence_masked_stop_gained(self):
        ''' test that get_mnv_consequence() works correctly
        '''
        
        var1 = make_vcf_line(extra='Codons=taT/taG').split('\t')
        var2 = make_vcf_line(extra='Codons=Tat/Cat').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        self.assertEqual(check_mnv_consequence(var1, var2, self.pattern),
            'masked_stop_gain_mnv')
    
    def test_check_mnv_consequence_alternate_residue(self):
        ''' test that get_mnv_consequence() works correctly
        '''
        
        var1 = make_vcf_line(extra='Codons=aaT/aaG').split('\t')
        var2 = make_vcf_line(extra='Codons=Aat/Cat').split('\t')
        
        var1 = parse_vcf_line(var1, self.Variant)
        var2 = parse_vcf_line(var2, self.Variant)
        
        self.assertEqual(check_mnv_consequence(var1, var2, self.pattern),
            'alternate_residue_mnv')
        
