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

from clinicalfilter.variant.info import Info
from clinicalfilter.variant.symbols import Symbols

class TestVariantInfoPy(unittest.TestCase):
    """  unit testing of the Info class
    """
    
    def setUp(self):
        """ define a default Info object
        """
        
        self.pops = ["AFR_AF", "AMR_AF", "ASN_AF", "DDD_AF", "EAS_AF",
            "ESP_AF", "EUR_AF", "MAX_AF", "SAS_AF", "UK10K_cohort_AF"]
        Info.set_populations(self.pops)
        
        # set up a SNV object, since SNV inherits Info
        self.info = Info("HGNC_ID=1001;CQ=missense_variant;random_tag")
    
    def tearDown(self):
        Info.set_populations([])
    
    def test_get_consequence(self):
        """ test that get_consequence works correctly
        """
        
        chrom, pos = '1', 1000
        info = Info('CQ=missense_variant;HGNC=TEST')
        alts = ('C',)
        
        # check that in the absence of any known conserved final exon positions,
        # the consequence is unchanged.
        self.assertEqual(info.get_consequences(chrom, pos, alts, []),
            [['missense_variant']])
        
        info = Info('CQ=missense_variant|stop_gained;HGNC=TEST|TEST2')
        self.assertEqual(info.get_consequences(chrom, pos, alts, []),
            [['missense_variant', 'stop_gained']])
        
    def test_get_consequence_last_base(self):
        '''check get_consequence() works with last base of exon changes
        '''
        
        chrom, pos = '1', 1000
        alts = ('C',)
        info = Info('CQ=missense_variant;HGNC=TEST')
        info.set_genes_and_consequence(chrom, pos, alts, [])
        
        # Now check that if the variant is at a position where it is a final
        # base in an exon with a conserved base, the consequence gets converted.
        info.last_base = set([("1", 1000)])
        self.assertEqual(info.get_consequences(chrom, pos, alts, []),
            [["conserved_exon_terminus_variant"]])
        
        # If we have a variant in multiple genes, check that it only alters the
        # missense/splice_region variants, and doesn't alter synonymous variants
        # (since these will be in transcripts where the variant is distant from
        # an exon boundary.)
        info = Info('CQ=missense_variant|synonymous_variant;HGNC=TEST|TEST1')
        info.set_genes_and_consequence(chrom, pos, alts, [])
        info.last_base = set([("1", 1000)])
        self.assertEqual(info.get_consequences(chrom, pos, alts, []),
            [["conserved_exon_terminus_variant", "synonymous_variant"]])
    
    def test_get_consequence_multiallelic(self):
        ''' test that get_consequence works correctly with multiple alleles
        '''
        
        chrom, pos = '1', 1000
        info = Info('CQ=missense_variant,synonymous_variant')
        alts = ('C', 'G')
        
        self.assertEqual(info.get_consequences(chrom, pos, alts, []),
            [['missense_variant'], ['synonymous_variant']])
        
    def test_get_consequence_multiallelic_with_masked(self):
        ''' test that get_consequence works correctly with multiple alleles
        '''
        
        chrom, pos = '1', 1000
        info = Info('CQ=missense_variant,synonymous_variant')
        alts = ('C', 'G')
        
        self.assertEqual(info.get_consequences(chrom, pos, alts, ['G']),
            [['missense_variant']])
    
    def test_parse_gene_symbols(self):
        """ test that parse_gene_symbols() works correctly
        """
        
        alts = ('C',)
        
        # check for when a HGNC key exists
        self.info["HGNC_ID"] = "A"
        genes = self.info.parse_gene_symbols(alts, [])
        self.assertEqual(genes, [Symbols(info={'HGNC_ID': 'A'}, idx=0)])
        
        # check for when a HGNC key doesn't exist
        del self.info["HGNC_ID"]
        genes = self.info.parse_gene_symbols(alts, [])
        self.assertEqual(genes, [Symbols(info={}, idx=0)])
        
        # check for multiple gene symbols
        self.info["HGNC_ID"] = "A|B|C"
        genes = self.info.parse_gene_symbols(alts, [])
        self.assertEqual(genes, [Symbols(info={'HGNC_ID': 'A|B|C'}, idx=0)])
        
        # check for multiple gene symbols, when some are missing
        self.info["HGNC_ID"] = "|.|C"
        genes = self.info.parse_gene_symbols(alts, [])
        self.assertEqual(genes, [Symbols(info={'HGNC_ID': '||C'}, idx=0)])
        
        # check for multiple gene symbols, when some missing symbols have
        # alternates in other symbol fields.
        self.info["HGNC_ID"] = ".|.|C"
        self.info["HGNC"] = "Z|.|C"
        genes = self.info.parse_gene_symbols(alts, [])
        self.assertEqual(genes, [Symbols(info={'HGNC_ID': '||C', 'HGNC': 'Z||C'}, idx=0)])
        
        # Check that including alternate symbols has the correct precendence
        # order. Note that doing this properly would require checking all of the
        # possible order combinations.
        self.info["HGNC_ID"] = ".|.|C"
        self.info["HGNC"] = "Z|.|C"
        self.info["SYMBOL"] = "A|.|C"
        genes = self.info.parse_gene_symbols(alts, [])
        self.assertEqual(genes, [Symbols(info={'HGNC_ID': '||C',
            'HGNC': 'Z||C', "SYMBOL": "A||C"}, idx=0)])
    
    def test_parse_gene_symbols_multi_alts(self):
        ''' check parse_gene_symbols() when we have multiple alleles
        '''
        
        info = Info('HGNC_ID=D,E;HGNC=D,E;SYMBOL=D,E;ENSG=D,E;ENST=D,E;ENSP=D,E;ENSR=D,E')
        alts = ('G', 'C')
        
        self.assertEqual(info.parse_gene_symbols(alts,  []),
            [Symbols(info={'HGNC_ID': 'D', 'HGNC': 'D', 'SYMBOL': 'D',
                'ENSG': 'D', 'ENST': 'D', 'ENSP': 'D', 'ENSR': 'D'}, idx=0),
            Symbols(info={'HGNC_ID': 'E', 'HGNC': 'E', 'SYMBOL': 'E',
                'ENSG': 'E', 'ENST': 'E', 'ENSP': 'E', 'ENSR': 'E'}, idx=0)])
        
        # if we have more alleles than the available symbols, we get an error
        # NOTE: this doesn't check if we have fewer alleles than symbols
        alts = ('G', 'T', 'C')
        with self.assertRaises(IndexError):
            self.info.parse_gene_symbols(alts, [])
        
    def test_parse_gene_symbols_multi_alts_multi_symbols(self):
        ''' check parse_gene_symbols() when we have multiple symbols per allele
        '''
        
        info = Info('HGNC_ID=D|X,E|Y;HGNC=D|X,E|Y;SYMBOL=D|X,E|Y;ENSG=D|X,E|Y;' \
            'ENST=D|X,E|Y;ENSP=D|X,E|Y;ENSR=D|X,E|Y')
        alts = ('G', 'C')
        
        self.assertEqual(info.parse_gene_symbols(alts, []),
            [Symbols(info={'HGNC_ID': 'D|X', 'HGNC': 'D|X', 'SYMBOL': 'D|X',
                'ENSG': 'D|X', 'ENST': 'D|X', 'ENSP': 'D|X', 'ENSR': 'D|X'}, idx=0),
            Symbols(info={'HGNC_ID': 'E|Y', 'HGNC': 'E|Y', 'SYMBOL': 'E|Y',
                'ENSG': 'E|Y', 'ENST': 'E|Y', 'ENSP': 'E|Y', 'ENSR': 'E|Y'}, idx=0)])
        
    def test_parse_gene_symbols_multi_alts_masked_alt(self):
        ''' check parse_gene_symbols() when we mask alt alleles
        '''
        
        info = Info('HGNC_ID=D|X,E|Y;HGNC=D|X,E|Y;SYMBOL=D|X,E|Y;ENSG=D|X,E|Y;' \
            'ENST=D|X,E|Y;ENSP=D|X,E|Y;ENSR=D|X,E|Y')
        alts = ('G', 'C')
        
        # mask one allele
        self.assertEqual(info.parse_gene_symbols(alts, ['C']),
            [Symbols(info={'HGNC_ID': 'D|X', 'HGNC': 'D|X', 'SYMBOL': 'D|X',
                'ENSG': 'D|X', 'ENST': 'D|X', 'ENSP': 'D|X', 'ENSR': 'D|X'}, idx=0)])
        
        # mask both alleles
        self.assertEqual(info.parse_gene_symbols(alts, ['C', 'G']),
            [])
    
    def test_parse_gene_symbols_missing_gene(self):
        ''' check the gene symbol is the genome pos when we lack any other info
        '''
        
        # remove the only possibly source of the gene symbol
        info = Info('')
        alts = ('C', )
        
        genes = info.parse_gene_symbols(alts, [])
        self.assertEqual(genes, [Symbols(info={}, idx=0)])
    
    def test_is_lof(self):
        """ test that is_lof() works correctly
        """
        
        # check that known LOF consensequence return True
        info = Info('CQ=stop_gained;HGNC=TEST')
        info.set_genes_and_consequence('1', 100, ('G'), [])
        self.assertTrue(info.is_lof())
        
        # check that known non-LOF consensequence returns False
        info = Info('CQ=missense_variant;HGNC=TEST')
        info.set_genes_and_consequence('1', 100, ('G'), [])
        self.assertFalse(info.is_lof())
        
        # check that null values return False
        info = Info('HGNC=TEST')
        info.set_genes_and_consequence('1', 100, ('G'), [])
        self.assertFalse(info.is_lof())
        
        # check when the variant overlaps multiple genes (so has multiple
        # gene symbols and consequences).
        info = Info('CQ=stop_gained|missense_variant;HGNC=ATRX|TTN')
        info.set_genes_and_consequence('1', 100, ('G'), [])
        
        self.assertTrue(info.is_lof())
        self.assertTrue(info.is_lof("ATRX"))
        self.assertFalse(info.is_lof("TTN"))
        
        # check that when we have a MNV, we can lose or gain a LOF annotation
        info.mnv_code = 'masked_stop_gain_mnv'
        self.assertFalse(info.is_lof("ATRX"))
        
        info.mnv_code = 'modified_stop_gained_mnv'
        self.assertTrue(info.is_lof("TTN"))
    
    def test_is_missense(self):
        """ test that is_missense() works correctly
        """
        
        # check that known missense equivalent consequence return True
        info = Info('CQ=missense_variant;HGNC=TEST')
        info.set_genes_and_consequence('1', 100, ('G'), [])
        self.assertTrue(info.is_missense(is_cnv=False))
        
        # check that known LoF equivalent consequence returns False
        info = Info('CQ=stop_gained;HGNC=TEST')
        info.set_genes_and_consequence('1', 100, ('G'), [])
        self.assertFalse(info.is_missense(is_cnv=False))
        
        # check that null values return False
        info = Info('HGNC=TEST')
        info.set_genes_and_consequence('1', 100, ('G'), [])
        self.assertFalse(info.is_missense(is_cnv=False))
        
        # check when the variant overlaps multiple genes (so has multiple
        # gene symbols and consequences).
        info = Info('CQ=missense_variant|synonymous_variant;HGNC=ATRX|TTN')
        info.set_genes_and_consequence('1', 100, ('G'), [])
        self.assertTrue(info.is_missense(is_cnv=False))
        self.assertTrue(info.is_missense(False, "ATRX"))
        self.assertFalse(info.is_missense(False, "TTN"))
        
        # check that when we have a MNV, we can lose or gain a LOF annotation
        info.mnv_code = 'modified_synonymous_mnv'
        self.assertFalse(info.is_missense(False, "ATRX"))
        
        info.mnv_code = 'modified_protein_altering_mnv'
        self.assertTrue(info.is_missense(False, "TTN"))
        
        # check that masked stop gained MNVs are converted to a missense
        info = Info('CQ=stop_gained;HGNC=ATRX')
        info.set_genes_and_consequence('1', 100, ('G'), [])
        info.mnv_code = 'masked_stop_gain_mnv'
        self.assertTrue(info.is_missense(False))
    
    def test_is_missense_cnv(self):
        ''' test that is_missense() works correctly for CNVs
        '''
        
        chrom, pos, alts, = '1', '15000000', ('G',)
        info = Info('HGNC=ATRX;CQ=coding_sequence_variant;random_tag')
        info.set_genes_and_consequence(chrom, pos, alts, [])
        
        self.assertTrue(info.is_missense(is_cnv=True))
        self.assertFalse(info.is_missense(is_cnv=False))
    
    def test_get_per_gene_consequence(self):
        """ test that get_per_gene_consequence works correctly
        """
        
        self.info.symbols = [Symbols(info={'HGNC': 'ATRX'}, idx=0)]
        self.info.consequence = [["missense_variant"]]
        
        self.assertEqual(self.info.get_per_gene_consequence(None), ["missense_variant"])
        self.assertEqual(self.info.get_per_gene_consequence("ATRX"), ["missense_variant"])
        self.assertEqual(self.info.get_per_gene_consequence("TEST"), [])
        
        # check a variant with consequences in multiple genes, that we only
        # pull out the consequencesquences for a single gene
        self.info.symbols = [Symbols(info={'HGNC': 'ATRX|TTN'}, idx=0)]
        self.info.consequence = [["missense_variant", "synonymous_variant"]]
        self.assertEqual(self.info.get_per_gene_consequence("ATRX"), ["missense_variant"])
        self.assertEqual(self.info.get_per_gene_consequence("TTN"), ["synonymous_variant"])
        
        # check a symbol where two symbols match, we only use the first consequence
        self.info.symbols = [Symbols(info={'HGNC': 'TEMP|ATRX|TEMP'}, idx=0)]
        self.info.consequence = [["splice_acceptor_variant", "missense_variant",
            "synonymous_variant"]]
        self.assertEqual(self.info.get_per_gene_consequence("TEMP"),
            ["splice_acceptor_variant"])
        
        # check a symbol with some None gene symbols
        self.info.symbols = [Symbols(info={'HGNC': '|ATRX|'}, idx=0)]
        self.info.consequence = [["splice_acceptor_variant", "missense_variant",
            "synonymous_variant"]]
        self.assertEqual(self.info.get_per_gene_consequence("ATRX"),
            ["missense_variant"])
    
    def test_get_allele_frequency(self):
        """ tests that number conversion works as expected
        """
        
        # single number returns that number
        self.assertEqual(self.info.get_allele_frequency("1"), 1)
        
        # two numbers return one number
        self.assertEqual(self.info.get_allele_frequency("1,1"), 1)
        
        # two numbers return the highest number
        self.assertEqual(self.info.get_allele_frequency("1,2"), 2)
        
        # number and string return the number
        self.assertEqual(self.info.get_allele_frequency("a,1"), 1)
        
        # single string value returns None
        self.assertEqual(self.info.get_allele_frequency("a"), None)
        
        # multiple string values return None
        self.assertEqual(self.info.get_allele_frequency("a,b"), None)
        
        # multiple string values return None
        self.assertEqual(self.info.get_allele_frequency(None), None)
    
    def test_is_number(self):
        """ tests that we can check if a value represents a number
        """
        
        self.assertEqual(self.info.is_number(None), False)
        self.assertEqual(self.info.is_number("5"), True)
        self.assertEqual(self.info.is_number("a"), False)
    
    def test_find_max_allele_frequency(self):
        """ test if the MAF finder operates correctly
        """
        
        # check for var without recorded MAF
        self.assertIsNone(self.info.find_max_allele_frequency())
        
        # check for single population
        self.info["MAX_AF"] = "0.005"
        self.assertEqual(self.info.find_max_allele_frequency(), 0.005)
        
        # check for two populations
        self.info["AFR_AF"] = "0.01"
        self.assertEqual(self.info.find_max_allele_frequency(), 0.01)
        
        # check for all populations
        pops = set(["AFR_AF", "AMR_AF", "ASN_AF", "DDD_AF", "EAS_AF", \
            "ESP_AF", "EUR_AF", "MAX_AF", "SAS_AF", "UK10K_cohort_AF"])
        for pop in pops:
            self.info[pop] = "0.05"
            self.assertEqual(self.info.find_max_allele_frequency(), 0.05)
        
        # make sure we can handle having None values
        self.info["AFR_AF"] = None
        self.assertEqual(self.info.find_max_allele_frequency(), 0.05)
    
    def test_find_max_allele_frequency_without_populations(self):
        ''' test if the MAF finder operates correctly when we haven't set any
        populations to check
        '''
        
        self.info["MAX_AF"] = "0.005"
        
        # this is a regression test for a problem that only occurs if the unit
        # tests are run in an order such that the populations might not have
        # been set in previous commits.
        Info.set_populations([])
        self.assertEqual(self.info.find_max_allele_frequency(), None)
        
        # reset the populations, so that other unit tests can also rely on the
        # populations being set
        Info.set_populations(self.pops)
