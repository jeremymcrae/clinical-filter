""" unit testing of the CNV class
"""

import unittest
from clinicalfilter.variant_cnv import CNV


class TestVariantCnvPy(unittest.TestCase):
    """
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
        
        tags = {"gene": ["HGNC", "VGN", "GN"], "consequence": \
            ["VCQ", "CQ"]}
        
        info = "HGNC=TEST;HGNC_ALL=TEST,OR5A1;CQ=missense_variant;CNSOLIDATE;WSCORE=0.5;CALLP=0.000;COMMONFORWARDS=0.000;MEANLR2=0.5;MADL2R=0.02;END=16000000;SVLEN=1000000"
        format_keys = "inheritance:DP"
        sample_values = "deNovo:50"
        
        self.var.add_info(info, tags)
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
        self.var.position = str(pseudoautosomal_region_start + 1000)
        self.var.info["END"] = str(pseudoautosomal_region_end - 1000)
        self.var.set_gender("F")
        
        self.var.alt_allele = "<DUP>"
        self.var.set_genotype()
        self.assertEqual(self.var.genotype, "DUP")
        self.assertEqual(self.var.get_inheritance_type(), "autosomal")
        
    def test_set_range(self):
        """ test that set_range() operates correctly
        """
        
        # check that range is set correctly under normal function
        self.var.position = "1000"
        self.var.info["END"] = "2000"
        self.var.set_range()
        self.assertEqual(self.var.range, ("1000", "2000"))
        
        # check that range is set correctly when no info available
        self.var.info = {}
        self.var.set_range()
        self.assertEqual(self.var.range, ("1000", "11000"))
    
    def test_fix_gene_IDs(self):
        """ test that fix_gene_IDs() works correctly
        """
        
        known_genes = {"TEST": {"start": "1000", "end": "2000", "chrom": "5"}}
        
        # make a CNV that will overlap with the known gene set
        self.var.gene = "TEST"
        self.var.position = "1000"
        self.var.info["END"] = "1500"
        
        # check that fixing gene names does not alter anything for a CNV in a 
        # single known gene
        self.var.fix_gene_IDs(known_genes)
        self.assertEqual(self.var.gene, "TEST")
        
        # check that fixing gene names does not alter names not in the gene dict
        self.var.gene = "TEST,TEST2"
        self.var.fix_gene_IDs(known_genes)
        self.assertEqual(self.var.gene, "TEST,TEST2")
        
        # check that fixing gene names drop name of genes where the name is in 
        # the known genes dict, and the CNV and gene do not overlap
        self.var.position = "900"
        self.var.info["END"] = "950"
        self.var.fix_gene_IDs(known_genes)
        self.assertEqual(self.var.gene, "TEST2")
        
        # check that when we do not have any known genes, the gene names are 
        # unaltered
        self.var.gene = "TEST,TEST2"
        self.var.fix_gene_IDs(None)
        self.assertEqual(self.var.gene, "TEST,TEST2")
    
    def test_add_gene_from_info_cnv(self):
        """ test that test_add_gene_from_info() works correctly
        """
        
        # check that HGNC_ALL takes precedence
        self.var.info["HGNC"] = "A"
        self.var.info["HGNC_ALL"] = "B"
        self.var.add_gene_from_info()
        self.assertEqual(self.var.gene, "B")
        
        # check that HGNC is used in the absence of HGNC_ALL
        del self.var.info["HGNC_ALL"]
        self.var.add_gene_from_info()
        self.assertEqual(self.var.gene, "A")
        
        # check that when HGNC and HGNC_ALL are undefined, we can still include
        # CNVs overlapping genes through NUMBERGENES > 0. 
        del self.var.info["HGNC"]
        
        # first test for NUMBERGENES = 0
        self.var.info["NUMBERGENES"] = 0
        self.var.add_gene_from_info()
        self.assertIsNone(self.var.gene)
        
        # and then make sure we are correct for NUMBERGENES > 0
        self.var.info["NUMBERGENES"] = 1
        self.var.add_gene_from_info()
        self.assertEqual(self.var.gene, ".")
        
        # finally check for no HGNC, HGNC_ALL, or NUMBERGENES
        del self.var.info["NUMBERGENES"]
        self.var.add_gene_from_info()
        self.assertIsNone(self.var.gene)
    
    def test_get_genes(self):
        """ test that get_genes() works correctly
        """
        
        self.var.gene = None
        self.assertEqual(self.var.get_genes(), [])
        
        self.var.gene = "TEST"
        self.assertEqual(self.var.get_genes(), ["TEST"])
        
        self.var.gene = "TEST1,TEST2"
        self.assertEqual(self.var.get_genes(), ["TEST1", "TEST2"])
        
        self.var.gene = "."
        self.assertEqual(self.var.get_genes(), ["."])
        
        self.var.gene = ","
        self.assertEqual(self.var.get_genes(), ["", ""])
    
    def test_fails_y_chrom_female(self):
        """ test that passes_filters() works correctly for female Y chrom CNVs
        """
        
        self.var.chrom = "Y"
        self.var.set_gender("F")
        
        filters = "temp"
        self.assertFalse(self.var.passes_filters(filters))


if __name__ == '__main__':
    unittest.main()


