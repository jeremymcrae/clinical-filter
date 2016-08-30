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
        known = {"ATRX": {"inheritance": {"Hemizygous": \
            {"Loss of function"}}, "start": "10000000", "chrom": "1", \
            "confirmed_status": {"Confirmed DD Gene"}, "end": "20000000"}}
        
        SNV.set_debug('1', 15000000)
        SNV.set_known_genes(known)
        
        # set up a SNV object, since SNV inherits Info
        info = "HGNC=ATRX;CQ=missense_variant;random_tag"
        self.var = SNV(chrom, pos, snp_id, ref, alt, filt, info=info)
    
    def test_set_consequence(self):
        """ test that set_consequence works correctly
        """
        
        # check that in the absence of any known conserved final exon positions,
        # the consequence is unchanged.
        self.var.set_consequence()
        self.assertEqual(self.var.consequence, ["missense_variant"])
        
        # Now check that if the variant is at a position where it is a final
        # base in an exon with a conserved base, the consequence gets converted.
        self.var.last_base = set([("1", 15000000)])
        self.var.set_consequence()
        self.assertEqual(self.var.consequence, ["conserved_exon_terminus_variant"])
        
        # If we have a variant in multiple genes, check that it only alters the
        # missense/splice_region variants, and doesn't alter synonymous variants
        # (since these will be in transcripts where the variant is distant from
        # an exon boundary.)
        self.var.info["CQ"] = "missense_variant|synonymous_variant"
        self.var.set_consequence()
        self.assertEqual(self.var.consequence, ["conserved_exon_terminus_variant", "synonymous_variant"])
    
    def test_set_gene_from_info(self):
        """ test that test_set_gene_from_info() works correctly
        """
        
        # check for when a HGNC key exists
        self.var.info["HGNC"] = "A"
        self.var.set_gene_from_info()
        self.assertEqual(self.var.genes, ["A"])
        
        # check for when a HGNC key doesn't exist
        del self.var.info["HGNC"]
        self.var.set_gene_from_info()
        self.assertIsNone(self.var.genes)
        
        # check for multiple gene symbols
        self.var.info["HGNC"] = "A|B|C"
        self.var.set_gene_from_info()
        self.assertEqual(self.var.genes, ["A", "B", "C"])
        
        # check for multiple gene symbols, when some are missing
        self.var.info["HGNC"] = "|.|C"
        self.var.set_gene_from_info()
        self.assertEqual(self.var.genes, [None, None, "C"])
        
        # check for multiple gene symbols, when some missing symbols have
        # alternates in other symbol fields.
        self.var.info["HGNC"] = ".|.|C"
        self.var.info["SYMBOL"] = "Z|.|C"
        self.var.set_gene_from_info()
        self.assertEqual(self.var.genes, ["Z", None, "C"])
        
        # Check that including alternate symbols has the correct precendence
        # order. Note that doing this properly would require checking all of the
        # possible order combinations.
        self.var.info["HGNC"] = ".|.|C"
        self.var.info["SYMBOL"] = "Z|.|C"
        self.var.info["ENSG"] = "A|.|C"
        self.var.set_gene_from_info()
        self.assertEqual(self.var.genes, ["Z", None, "C"])
    
    def test_set_gene_from_info_missing_gene(self):
        ''' check the gene symbol is the genome pos when we lack any other info
        '''
        
        # remove the known genes, an=y previously set gene info
        self.var.known_genes = None
        self.var.genes = None
        
        # remove the only possibly source of the gene symbol
        del self.var.info["HGNC"]
        
        self.var.set_gene_from_info()
        self.assertEqual(self.var.genes, ["1:15000000"])
    
    def test_is_lof(self):
        """ test that is_lof() works correctly
        """
        
        # check that known LOF consensequence return True
        self.var.consequence = ["stop_gained"]
        self.assertTrue(self.var.is_lof())
        
        # check that known non-LOF consensequence returns False
        self.var.consequence = ["missense_variant"]
        self.assertFalse(self.var.is_lof())
        
        # check that null values return False
        self.var.consequence = None
        self.assertFalse(self.var.is_lof())
        
        # check when the variant overlaps multiple genes (so has multiple
        # gene symbols and consequences).
        self.var.consequence = ["stop_gained", "missense_variant"]
        self.var.genes = ["ATRX", "TTN"]
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
        self.var.consequence = ["missense_variant"]
        self.assertTrue(self.var.is_missense())
        
        # check that known LoF equivalent consequence returns False
        self.var.consequence = ["stop_gained"]
        self.assertFalse(self.var.is_missense())
        
        # check that null values return False
        self.var.consequence = None
        self.assertFalse(self.var.is_missense())
        
        # check when the variant overlaps multiple genes (so has multiple
        # gene symbols and consequences).
        self.var.consequence = ["missense_variant", "synonymous_variant"]
        self.var.genes = ["ATRX", "TTN"]
        self.assertTrue(self.var.is_missense())
        self.assertTrue(self.var.is_missense("ATRX"))
        self.assertFalse(self.var.is_missense("TTN"))
        
        # check that when we have a MNV, we can lose or gain a LOF annotation
        self.var.mnv_code = 'modified_synonymous_mnv'
        self.assertFalse(self.var.is_missense("ATRX"))
        
        self.var.mnv_code = 'modified_protein_altering_mnv'
        self.assertTrue(self.var.is_missense("TTN"))
    
    def test_correct_multiple_alt(self):
        """ test that correct_multiple_alt works correctly
        """
        
        # define the number of alleles and consequences for multiple alleles
        self.var.info["AC"] = "1,1"
        cq = ["missense_variant,splice_acceptor_variant"]
        
        # check with alts that fall in one gene
        self.var.info["HGNC"] = "ATRX,ATRX"
        self.var.set_gene_from_info()
        self.assertEqual(self.var.correct_multiple_alt(cq),
            (['splice_acceptor_variant'], ['ATRX'], None))
        
        # check with alts that fall in multiple genes
        cq = ["missense_variant|regulatory_region_variant,stop_gained|splice_acceptor_variant"]
        self.var.info["HGNC"] = "ATRX|TTN,ATRX|TTN"
        self.var.set_gene_from_info()
        self.assertEqual(self.var.correct_multiple_alt(cq),
            (['stop_gained', 'splice_acceptor_variant'], ['ATRX', 'TTN'], None))
        
        # check a cq that has already been split by "|" (ie by gene)
        cq = ["missense_variant", "regulatory_region_variant,stop_gained",
            "splice_acceptor_variant"]
        self.var.set_gene_from_info()
        self.assertEqual(self.var.correct_multiple_alt(cq),
            (['stop_gained', 'splice_acceptor_variant'], ['ATRX', 'TTN'], None))
        
        # check that if the proband has a zero count for an allele, then we
        # disregard the consequences and HGNC symbols for that allele
        self.var.info["AC"] = "1,0"
        self.var.set_gene_from_info()
        self.assertEqual(self.var.correct_multiple_alt(cq),
            (['missense_variant', 'regulatory_region_variant'], ['ATRX', 'TTN'], None))
        
        # revert the allele counts, but drop the HGNC symbol, and make sure the
        # HGNC symbol returned is None
        self.var.info["AC"] = "1,1"
        del self.var.info["HGNC"]
        self.var.set_gene_from_info()
        self.assertEqual(self.var.correct_multiple_alt(cq),
            (['stop_gained', 'splice_acceptor_variant'], [], None))
    
    def test_correct_multiple_alt_enst(self):
        ''' check the impact of correcting multiple alt cq on ENST values
        '''
        
        # define the number of alleles and consequences for multiple alleles
        self.var.info["AC"] = "1,1"
        cq = ["missense_variant,splice_acceptor_variant"]
        
        # check with alts that fall in one gene
        self.var.info["HGNC"] = "ATRX,ATRX"
        self.var.info["ENST"] = "ENST100,ENST200"
        self.var.set_gene_from_info()
        
        self.assertEqual(self.var.correct_multiple_alt(cq),
            (['splice_acceptor_variant'], ['ATRX'], 'ENST100'))
    
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
        
        self.var.genes = ["ATRX"]
        self.var.consequence = ["missense_variant"]
        
        self.assertEqual(self.var.get_per_gene_consequence(None), ["missense_variant"])
        self.assertEqual(self.var.get_per_gene_consequence("ATRX"), ["missense_variant"])
        self.assertEqual(self.var.get_per_gene_consequence("TEST"), [])
        
        # check a variant with consequences in multiple genes, that we only
        # pull out the consequencesquences for a single gene
        self.var.genes = ["ATRX", "TTN"]
        self.var.consequence = ["missense_variant", "synonymous_variant"]
        self.assertEqual(self.var.get_per_gene_consequence("ATRX"), ["missense_variant"])
        self.assertEqual(self.var.get_per_gene_consequence("TTN"), ["synonymous_variant"])
        
        # check a symbol where two symbols match
        self.var.genes = ["TEMP", "ATRX", "TEMP"]
        self.var.consequence = ["splice_acceptor_variant", "missense_variant", \
            "synonymous_variant"]
        self.assertEqual(self.var.get_per_gene_consequence("TEMP"), \
            ["splice_acceptor_variant", "synonymous_variant"])
        
        # check a symbol with some None gene symbols
        self.var.genes = [None, "ATRX", None]
        self.var.consequence = ["splice_acceptor_variant", "missense_variant", \
            "synonymous_variant"]
        self.assertEqual(self.var.get_per_gene_consequence("ATRX"), \
            ["missense_variant"])
        
        # check that the earlier VCFs with single consequences but multiple
        # symbols from HGNC_ALL give the same consequence for all genes.
        info = "HGNC_ALL=ATRX&TTN;CQ=missense_variant;random_tag"
        del self.var.info["HGNC"]
        self.var.genes = None
        self.var.add_info(info)
        
        self.assertEqual(self.var.get_per_gene_consequence("ATRX"), \
            ["missense_variant"])
        self.assertEqual(self.var.get_per_gene_consequence("TTN"), \
            ["missense_variant"])
        
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
