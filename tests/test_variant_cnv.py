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
from clinicalfilter.variant.cnv import CNV


class TestVariantCnvPy(unittest.TestCase):
    """ unit testing of the CNV class
    """
    
    def setUp(self):
        """ define a default VcfInfo object
        """
        
        chrom = "1"
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "<DUP>"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        self.var = CNV(chrom, pos, snp_id, ref, alt, filt)
        
        info = "HGNC=TEST;HGNC_ALL=TEST,OR5A1;CQ=missense_variant;CNSOLIDATE;WSCORE=0.5;CALLP=0.000;COMMONFORWARDS=0.000;MEANLR2=0.5;MADL2R=0.02;END=16000000;SVLEN=1000000"
        format_keys = "inheritance:DP"
        sample_values = "deNovo:50"
        
        self.var.add_info(info)
        self.var.add_format(format_keys, sample_values)
        self.var.set_gender("F")
    
    def test_set_genotype(self):
        """ test that set_genotype() operates correctly
        """
        
        # check that DUPs are set correctly
        self.var.alt_allele = "<DUP>"
        self.var.set_genotype()
        self.assertEqual(self.var.genotype, "DUP")
        
        # check that DELs are set correctly
        self.var.alt_allele = "<DEL>"
        self.var.set_genotype()
        self.assertEqual(self.var.genotype, "DEL")
        
        # check that other genotypes raise an error
        self.var.alt_allele = "G"
        with self.assertRaises(ValueError):
            self.var.set_genotype()
        
        # and check that we raise an error for female Y chrom CNVs
        self.var.chrom = "Y"
        self.var.set_gender("F")
        with self.assertRaises(ValueError):
            self.var.set_genotype()
    
    def test_set_genotype_pseudoautosomal(self):
        """ test that set_genotype() works correctly in pseudoautosomal regions
        """
        
        pseudoautosomal_region_start = 60002
        pseudoautosomal_region_end = 2699520
        
        # set a CNV that lies within a pseudoautosomal region
        self.var.chrom = "X"
        self.var.position = pseudoautosomal_region_start + 1000
        self.var.info["END"] = pseudoautosomal_region_end - 1000
        self.var.set_gender("F")
        
        self.var.alt_allele = "<DUP>"
        self.var.set_genotype()
        self.assertEqual(self.var.genotype, "DUP")
        self.assertEqual(self.var.get_inheritance_type(), "autosomal")
        
    def test_get_range(self):
        """ test that get_range() operates correctly
        """
        
        # check that range is set correctly under normal function
        self.var.position = 1000
        self.var.info["END"] = "2000"
        self.assertEqual(self.var.get_range(), (1000, 2000))
        
        # check that range is set correctly when no info available
        self.var.info = {}
        self.assertEqual(self.var.get_range(), (1000, 11000))
    
    def test_fix_gene_IDs(self):
        """ test that fix_gene_IDs() works correctly
        """
        
        self.var.known_genes = {"TEST": {"start": 1000, "end": 2000, "chrom": "5"}}
        
        # make a CNV that will overlap with the known gene set
        self.var.genes = ["TEST"]
        self.var.position = 1000
        self.var.info["END"] = "1500"
        
        # check that fixing gene names does not alter anything for a CNV in a
        # single known gene
        self.var.fix_gene_IDs()
        self.assertEqual(self.var.genes, ["TEST"])
        
        # check that fixing gene names does not alter names not in the gene dict
        self.var.genes = ["TEST", "TEST2"]
        self.var.fix_gene_IDs()
        self.assertEqual(self.var.genes, ["TEST", "TEST2"])
        
        # check that fixing gene names drop name of genes where the name is in
        # the known genes dict, and the CNV and gene do not overlap
        self.var.position = 900
        self.var.info["END"] = "950"
        self.var.fix_gene_IDs()
        self.assertEqual(self.var.genes, [".", "TEST2"])
        
        # check that when we do not have any known genes, the gene names are
        # unaltered
        self.var.genes = ["TEST", "TEST2"]
        self.var.known_genes = None
        self.var.fix_gene_IDs()
        self.assertEqual(self.var.genes, ["TEST", "TEST2"])
    
    def test_set_gene_from_info_cnv(self):
        """ test that set_add_gene_from_info() works correctly
        """
        
        # make sure the known genes are None, otherwise sometimes the values
        # from test_variant_info.py unit tests can bleed through. I'm not sure
        # why!
        self.var.known_genes = None
        
        # check that HGNC takes precedence
        self.var.info["HGNC"] = "A"
        self.var.info["HGNC_ALL"] = "B"
        self.var.set_gene_from_info()
        self.assertEqual(self.var.genes, ["A"])
        
        # check that HGNC is used in the absence of HGNC_ALL
        del self.var.info["HGNC"]
        self.var.set_gene_from_info()
        self.assertEqual(self.var.genes, ["B"])
        
        # check that when HGNC and HGNC_ALL are undefined, we can still include
        # CNVs overlapping genes through NUMBERGENES > 0.
        del self.var.info["HGNC_ALL"]
        
        # first test for NUMBERGENES = 0
        self.var.info["NUMBERGENES"] = 0
        self.var.set_gene_from_info()
        self.assertIsNone(self.var.genes)
        
        # and then make sure we are correct for NUMBERGENES > 0
        self.var.info["NUMBERGENES"] = 1
        self.var.set_gene_from_info()
        self.assertEqual(self.var.genes, ["."])
        
        # finally check for no HGNC, HGNC_ALL, or NUMBERGENES
        del self.var.info["NUMBERGENES"]
        self.var.set_gene_from_info()
        self.assertEqual(self.var.genes, ["1:15000000"])
    
    def test_get_genes(self):
        """ test that get_genes() works correctly
        """
        
        self.var.genes = None
        self.assertEqual(self.var.get_genes(), [])
        
        self.var.genes = ["TEST"]
        self.assertEqual(self.var.get_genes(), ["TEST"])
        
        self.var.genes = ["TEST1", "TEST2"]
        self.assertEqual(self.var.get_genes(), ["TEST1", "TEST2"])
        
        self.var.genes = ["."]
        self.assertEqual(self.var.get_genes(), ["."])
    
    def test_fails_y_chrom_female(self):
        """ test that passes_filters() works correctly for female Y chrom CNVs
        """
        
        self.var.chrom = "Y"
        self.var.set_gender("F")
        
        self.assertFalse(self.var.passes_filters())


if __name__ == '__main__':
    unittest.main()
