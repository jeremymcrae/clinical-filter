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

from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.variant.symbols import Symbols

class TestVariantInfoPy(unittest.TestCase):
    """  unit testing of the VcfInfo class
    """
    
    def setUp(self):
        """ define a default VcfInfo object
        """
        
        chrom = "1"
        pos = "15000000"
        snp_id = "CM00001"
        ref = "A"
        alt = "G"
        filt = "PASS"
        
        # here are the default filtering criteria, as loaded into python
        known = {"1001": {"inh": {'symbol': 'ATRX', "Hemizygous":
            {"Loss of function"}}, "start": "10000000", "chrom": "1",
            "confirmed_status": {"confirmed dd gene"}, "end": "20000000"}}
        
        self.pops = ["AFR_AF", "AMR_AF", "ASN_AF", "DDD_AF", "EAS_AF",
            "ESP_AF", "EUR_AF", "MAX_AF", "SAS_AF", "UK10K_cohort_AF"]
        
        SNV.set_debug('1', 15000000)
        SNV.set_known_genes(known)
        SNV.set_populations(self.pops)
        
        # set up a SNV object, since SNV inherits Info
        info = "HGNC_ID=1001;CQ=missense_variant;random_tag"
        self.var = SNV(chrom, pos, snp_id, ref, alt, filt, info=info)
    
    def test_get_consequence(self):
        """ test that get_consequence works correctly
        """
        
        info = {'CQ': 'missense_variant'}
        alt_alleles = ('C',)
        
        # check that in the absence of any known conserved final exon positions,
        # the consequence is unchanged.
        self.assertEqual(self.var.get_consequences(info, alt_alleles, []),
            [['missense_variant']])
        
        info = {'CQ': 'missense_variant|stop_gained'}
        self.assertEqual(self.var.get_consequences(info, alt_alleles, []),
            [['missense_variant', 'stop_gained']])
        
    def test_get_consequence_last_base(self):
        '''check get_consequence() works with last base of exon changes
        '''
        
        info = {'CQ': 'missense_variant'}
        alt_alleles = ('C',)
        
        # Now check that if the variant is at a position where it is a final
        # base in an exon with a conserved base, the consequence gets converted.
        self.var.last_base = set([("1", 15000000)])
        self.assertEqual(self.var.get_consequences(info, alt_alleles, []),
            [["conserved_exon_terminus_variant"]])
        
        # If we have a variant in multiple genes, check that it only alters the
        # missense/splice_region variants, and doesn't alter synonymous variants
        # (since these will be in transcripts where the variant is distant from
        # an exon boundary.)
        info["CQ"] = "missense_variant|synonymous_variant"
        self.assertEqual(self.var.get_consequences(info, alt_alleles, []),
            [["conserved_exon_terminus_variant", "synonymous_variant"]])
    
    def test_get_consequence_multiallelic(self):
        ''' test that get_consequence works correctly with multiple alleles
        '''
        
        info = {'CQ': 'missense_variant,synonymous_variant'}
        alt_alleles = ('C', 'G')
        
        self.assertEqual(self.var.get_consequences(info, alt_alleles, []),
            [['missense_variant'], ['synonymous_variant']])
        
    def test_get_consequence_multiallelic_with_masked(self):
        ''' test that get_consequence works correctly with multiple alleles
        '''
        
        info = {'CQ': 'missense_variant,synonymous_variant'}
        alt_alleles = ('C', 'G')
        
        self.assertEqual(self.var.get_consequences(info, alt_alleles, ['G']),
            [['missense_variant']])
    
    def test_get_gene_from_info(self):
        """ test that test_get_gene_from_info() works correctly
        """
        
        # check for when a HGNC key exists
        self.var.info["HGNC_ID"] = "A"
        genes = self.var.get_gene_from_info(self.var.info, self.var.alt_alleles, [])
        self.assertEqual(genes, [Symbols(info={'HGNC_ID': 'A'}, idx=0)])
        
        # check for when a HGNC key doesn't exist
        del self.var.info["HGNC_ID"]
        genes = self.var.get_gene_from_info(self.var.info, self.var.alt_alleles, [])
        self.assertEqual(genes, [Symbols(info={}, idx=0)])
        
        # check for multiple gene symbols
        self.var.info["HGNC_ID"] = "A|B|C"
        genes = self.var.get_gene_from_info(self.var.info, self.var.alt_alleles, [])
        self.assertEqual(genes, [Symbols(info={'HGNC_ID': 'A|B|C'}, idx=0)])
        
        # check for multiple gene symbols, when some are missing
        self.var.info["HGNC_ID"] = "|.|C"
        genes = self.var.get_gene_from_info(self.var.info, self.var.alt_alleles, [])
        self.assertEqual(genes, [Symbols(info={'HGNC_ID': '||C'}, idx=0)])
        
        # check for multiple gene symbols, when some missing symbols have
        # alternates in other symbol fields.
        self.var.info["HGNC_ID"] = ".|.|C"
        self.var.info["HGNC"] = "Z|.|C"
        genes = self.var.get_gene_from_info(self.var.info, self.var.alt_alleles, [])
        self.assertEqual(genes, [Symbols(info={'HGNC_ID': '||C', 'HGNC': 'Z||C'}, idx=0)])
        
        # Check that including alternate symbols has the correct precendence
        # order. Note that doing this properly would require checking all of the
        # possible order combinations.
        self.var.info["HGNC_ID"] = ".|.|C"
        self.var.info["HGNC"] = "Z|.|C"
        self.var.info["SYMBOL"] = "A|.|C"
        genes = self.var.get_gene_from_info(self.var.info, self.var.alt_alleles, [])
        self.assertEqual(genes, [Symbols(info={'HGNC_ID': '||C',
            'HGNC': 'Z||C', "SYMBOL": "A||C"}, idx=0)])
    
    def test_get_gene_from_info_multi_alts(self):
        ''' check check get_gene_from_info() when we have multiple alleles
        '''
        
        info = self.var.info
        alt_alleles = ('G', 'C')
        info['HGNC_ID'] = 'D,E'
        info['HGNC'] = 'D,E'
        info['SYMBOL'] = 'D,E'
        info['ENSG'] = 'D,E'
        info['ENST'] = 'D,E'
        info['ENSP'] = 'D,E'
        info['ENSR'] = 'D,E'
        
        self.assertEqual(self.var.get_gene_from_info(info, alt_alleles,  []),
            [Symbols(info={'HGNC_ID': 'D', 'HGNC': 'D', 'SYMBOL': 'D',
                'ENSG': 'D', 'ENST': 'D', 'ENSP': 'D', 'ENSR': 'D'}, idx=0),
            Symbols(info={'HGNC_ID': 'E', 'HGNC': 'E', 'SYMBOL': 'E',
                'ENSG': 'E', 'ENST': 'E', 'ENSP': 'E', 'ENSR': 'E'}, idx=0)])
        
        # if we have more alleles than the available symbols, we get an error
        # NOTE: this doesn't check if we have fewer alleles than symbols
        alt_alleles = ('G', 'T', 'C')
        with self.assertRaises(IndexError):
            self.var.get_gene_from_info(info, alt_alleles, [])
        
    def test_get_gene_from_info_multi_alts_multi_symbols(self):
        ''' check get_gene_from_info() when we have multiple symbols per allele
        '''
        
        info = self.var.info
        alt_alleles = ('G', 'C')
        info['HGNC_ID'] = 'D|X,E|Y'
        info['HGNC'] = 'D|X,E|Y'
        info['SYMBOL'] = 'D|X,E|Y'
        info['ENSG'] = 'D|X,E|Y'
        info['ENST'] = 'D|X,E|Y'
        info['ENSP'] = 'D|X,E|Y'
        info['ENSR'] = 'D|X,E|Y'
        
        self.assertEqual(self.var.get_gene_from_info(info, alt_alleles, []),
            [Symbols(info={'HGNC_ID': 'D|X', 'HGNC': 'D|X', 'SYMBOL': 'D|X',
                'ENSG': 'D|X', 'ENST': 'D|X', 'ENSP': 'D|X', 'ENSR': 'D|X'}, idx=0),
            Symbols(info={'HGNC_ID': 'E|Y', 'HGNC': 'E|Y', 'SYMBOL': 'E|Y',
                'ENSG': 'E|Y', 'ENST': 'E|Y', 'ENSP': 'E|Y', 'ENSR': 'E|Y'}, idx=0)])
        
    def test_get_gene_from_info_multi_alts_masked_alt(self):
        ''' check get_gene_from_info() when we mask alt alleles
        '''
        
        info = self.var.info
        alt_alleles = ('G', 'C')
        info['HGNC_ID'] = 'D|X,E|Y'
        info['HGNC'] = 'D|X,E|Y'
        info['SYMBOL'] = 'D|X,E|Y'
        info['ENSG'] = 'D|X,E|Y'
        info['ENST'] = 'D|X,E|Y'
        info['ENSP'] = 'D|X,E|Y'
        info['ENSR'] = 'D|X,E|Y'
        
        # mask one allele
        self.assertEqual(self.var.get_gene_from_info(info, alt_alleles, ['C']),
            [Symbols(info={'HGNC_ID': 'D|X', 'HGNC': 'D|X', 'SYMBOL': 'D|X',
                'ENSG': 'D|X', 'ENST': 'D|X', 'ENSP': 'D|X', 'ENSR': 'D|X'}, idx=0)])
        
        # mask both alleles
        self.assertEqual(self.var.get_gene_from_info(info, alt_alleles, ['C', 'G']),
            [])
    
    def test_get_gene_from_info_missing_gene(self):
        ''' check the gene symbol is the genome pos when we lack any other info
        '''
        
        # remove the known genes, any previously set gene info
        self.var.known_genes = None
        self.var.genes = None
        
        # remove the only possibly source of the gene symbol
        del self.var.info["HGNC_ID"]
        
        genes = self.var.get_gene_from_info(self.var.info, self.var.alt_alleles, [])
        self.assertEqual(genes, [Symbols(info={}, idx=0)])
    
    def test_is_lof(self):
        """ test that is_lof() works correctly
        """
        
        # check that known LOF consensequence return True
        self.var.consequence = [["stop_gained"]]
        self.assertTrue(self.var.is_lof())
        
        # check that known non-LOF consensequence returns False
        self.var.consequence = [["missense_variant"]]
        self.assertFalse(self.var.is_lof())
        
        # check that null values return False
        self.var.consequence = None
        self.assertFalse(self.var.is_lof())
        
        # check when the variant overlaps multiple genes (so has multiple
        # gene symbols and consequences).
        self.var.consequence = [["stop_gained", "missense_variant"]]
        
        self.var.genes = [Symbols(info={'HGNC': 'ATRX|TTN'}, idx=0)]
        self.assertTrue(self.var.is_lof())
        self.assertTrue(self.var.is_lof("ATRX"))
        self.assertFalse(self.var.is_lof("TTN"))
        
        # check that when we have a MNV, we can lose or gain a LOF annotation
        self.var.mnv_code = 'masked_stop_gain_mnv'
        self.assertFalse(self.var.is_lof("ATRX"))
        
        self.var.mnv_code = 'modified_stop_gained_mnv'
        self.assertTrue(self.var.is_lof("TTN"))
    
    def test_is_missense(self):
        """ test that is_missense() works correctly
        """
        
        # check that known missense equivalent consequence return True
        self.var.consequence = [["missense_variant"]]
        self.assertTrue(self.var.is_missense())
        
        # check that known LoF equivalent consequence returns False
        self.var.consequence = [["stop_gained"]]
        self.assertFalse(self.var.is_missense())
        
        # check that null values return False
        self.var.consequence = None
        self.assertFalse(self.var.is_missense())
        
        # check when the variant overlaps multiple genes (so has multiple
        # gene symbols and consequences).
        self.var.consequence = [["missense_variant", "synonymous_variant"]]
        self.var.genes = [Symbols(info={'HGNC': 'ATRX|TTN'}, idx=0)]
        self.assertTrue(self.var.is_missense())
        self.assertTrue(self.var.is_missense("ATRX"))
        self.assertFalse(self.var.is_missense("TTN"))
        
        # check that when we have a MNV, we can lose or gain a LOF annotation
        self.var.mnv_code = 'modified_synonymous_mnv'
        self.assertFalse(self.var.is_missense("ATRX"))
        
        self.var.mnv_code = 'modified_protein_altering_mnv'
        self.assertTrue(self.var.is_missense("TTN"))
        
        # check that masked stop gained MNVs are converted to a missense
        self.var.consequence = [["stop_gained"]]
        self.var.mnv_code = 'masked_stop_gain_mnv'
        self.assertTrue(self.var.is_missense())
    
    def test_is_missense_cnv(self):
        ''' test that is_missense() works correctly for CNVs
        '''
        
        chrom, pos, var_id = '1', '15000000', '.'
        ref, alt, status = 'A', 'G', 'PASS'
        info = 'HGNC=ATRX;CQ=coding_sequence_variant;random_tag'
        
        snv = SNV(chrom, pos, var_id, ref, alt, status, info=info)
        cnv = CNV(chrom, pos, var_id, ref, alt, status, info=info)
        
        self.assertTrue(cnv.is_missense())
        self.assertFalse(snv.is_missense())
    
    def test_get_most_severe_consequence(self):
        """ test that get_most_severe_consequence works correctly
        """
        
        # check for the most simple list
        cq = ["missense_variant", "splice_acceptor_variant"]
        self.assertEqual(self.var.get_most_severe_consequence(cq), \
            "splice_acceptor_variant")
        
        # check for a single-entry list
        cq = ["missense_variant"]
        self.assertEqual(self.var.get_most_severe_consequence(cq), "missense_variant")
        
        # check for lists of lists per allele
        cq_per_allele = [["synonymous_variant", "splice_donor_variant"], \
            ["missense_variant", "regulatory_region_variant"]]
        self.assertEqual(self.var.get_most_severe_consequence(cq_per_allele), \
            ["missense_variant", "splice_donor_variant"])
    
    def test_get_per_gene_consequence(self):
        """ test that get_per_gene_consequence works correctly
        """
        
        self.var.genes = [Symbols(info={'HGNC': 'ATRX'}, idx=0)]
        self.var.consequence = [["missense_variant"]]
        
        self.assertEqual(self.var.get_per_gene_consequence(None), ["missense_variant"])
        self.assertEqual(self.var.get_per_gene_consequence("ATRX"), ["missense_variant"])
        self.assertEqual(self.var.get_per_gene_consequence("TEST"), [])
        
        # check a variant with consequences in multiple genes, that we only
        # pull out the consequencesquences for a single gene
        self.var.genes = [Symbols(info={'HGNC': 'ATRX|TTN'}, idx=0)]
        self.var.consequence = [["missense_variant", "synonymous_variant"]]
        self.assertEqual(self.var.get_per_gene_consequence("ATRX"), ["missense_variant"])
        self.assertEqual(self.var.get_per_gene_consequence("TTN"), ["synonymous_variant"])
        
        # check a symbol where two symbols match, we only use the first consequence
        self.var.genes = [Symbols(info={'HGNC': 'TEMP|ATRX|TEMP'}, idx=0)]
        self.var.consequence = [["splice_acceptor_variant", "missense_variant",
            "synonymous_variant"]]
        self.assertEqual(self.var.get_per_gene_consequence("TEMP"),
            ["splice_acceptor_variant"])
        
        # check a symbol with some None gene symbols
        self.var.genes = [Symbols(info={'HGNC': '|ATRX|'}, idx=0)]
        self.var.consequence = [["splice_acceptor_variant", "missense_variant",
            "synonymous_variant"]]
        self.assertEqual(self.var.get_per_gene_consequence("ATRX"),
            ["missense_variant"])
    
    def test_get_low_depth_alleles(self):
        ''' test that get_low_depth_alleles() works correctly
        '''
        
        # check with a single allele whre it is non-zero
        info = {'AC': '10'}
        alt_alleles = ('C', )
        self.assertEqual(self.var.get_low_depth_alleles(info, alt_alleles), [])
        
        # check with a single allele with zero count
        info = {'AC': '0'}
        alt_alleles = ('C', )
        self.assertEqual(self.var.get_low_depth_alleles(info, alt_alleles), ['C'])
        
        # check with multiallelic, where both are nonzero
        info = {'AC': '10,10'}
        alt_alleles = ('C', 'G')
        self.assertEqual(self.var.get_low_depth_alleles(info, alt_alleles), [])
        
        # check with multiallelic, where one a has zero count
        info = {'AC': '10,0'}
        alt_alleles = ('C', 'G')
        self.assertEqual(self.var.get_low_depth_alleles(info, alt_alleles), ['G'])
    
    def test_get_low_depth_alleles_indels(self):
        ''' test that get_low_depth_alleles() works correctly for indels
        '''
        # check this doesn't affect SNVs
        info = {'AC': '10,1'}
        alt_alleles = ('C', 'G')
        self.assertEqual(self.var.get_low_depth_alleles(info, alt_alleles), [])
        
        # check it works when the ref allele is > 1 bp
        self.var.ref_allele = 'GG'
        self.assertEqual(self.var.get_low_depth_alleles(info, alt_alleles), ['G'])
        
        # check it doesn't filter at 2+ reads
        info = {'AC': '10,2'}
        self.assertEqual(self.var.get_low_depth_alleles(info, alt_alleles), [])
        
        # check it catches when alt allele indicates an indel
        self.var.ref_allele = 'G'
        info = {'AC': '10,1'}
        alt_alleles = ('C', 'GG')
        self.assertEqual(self.var.get_low_depth_alleles(info, alt_alleles), ['GG'])
        
    def test_get_allele_frequency(self):
        """ tests that number conversion works as expected
        """
        
        # single number returns that number
        self.assertEqual(self.var.get_allele_frequency("1"), 1)
        
        # two numbers return one number
        self.assertEqual(self.var.get_allele_frequency("1,1"), 1)
        
        # two numbers return the highest number
        self.assertEqual(self.var.get_allele_frequency("1,2"), 2)
        
        # number and string return the number
        self.assertEqual(self.var.get_allele_frequency("a,1"), 1)
        
        # single string value returns None
        self.assertEqual(self.var.get_allele_frequency("a"), None)
        
        # multiple string values return None
        self.assertEqual(self.var.get_allele_frequency("a,b"), None)
        
        # multiple string values return None
        self.assertEqual(self.var.get_allele_frequency(None), None)
    
    def test_is_number(self):
        """ tests that we can check if a value represents a number
        """
        
        self.assertEqual(self.var.is_number(None), False)
        self.assertEqual(self.var.is_number("5"), True)
        self.assertEqual(self.var.is_number("a"), False)
    
    def test_find_max_allele_frequency(self):
        """ test if the MAF finder operates correctly
        """
        
        # check for var without recorded MAF
        self.assertIsNone(self.var.find_max_allele_frequency())
        
        # check for single population
        self.var.info["MAX_AF"] = "0.005"
        self.assertEqual(self.var.find_max_allele_frequency(), 0.005)
        
        # check for two populations
        self.var.info["AFR_AF"] = "0.01"
        self.assertEqual(self.var.find_max_allele_frequency(), 0.01)
        
        # check for all populations
        pops = set(["AFR_AF", "AMR_AF", "ASN_AF", "DDD_AF", "EAS_AF", \
            "ESP_AF", "EUR_AF", "MAX_AF", "SAS_AF", "UK10K_cohort_AF"])
        for pop in pops:
            self.var.info[pop] = "0.05"
            self.assertEqual(self.var.find_max_allele_frequency(), 0.05)
        
        # make sure we can handle having None values
        self.var.info["AFR_AF"] = None
        self.assertEqual(self.var.find_max_allele_frequency(), 0.05)
    
    def test_find_max_allele_frequency_without_populations(self):
        ''' test if the MAF finder operates correctly when we haven't set any
        populations to check
        '''
        
        self.var.info["MAX_AF"] = "0.005"
        
        # this is a regression test for a problem that only occurs if the unit
        # tests are run in an order such that the populations might not have
        # been set in previous commits.
        SNV.set_populations([])
        self.assertEqual(self.var.find_max_allele_frequency(), None)
        
        # reset the populations, so that other unit tests can also rely on the
        # populations being set
        SNV.set_populations(self.pops)
