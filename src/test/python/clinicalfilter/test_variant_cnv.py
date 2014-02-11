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
        qual = "50"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        self.var = CNV(chrom, pos, snp_id, ref, alt, qual, filt)
        
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
        del self.var.info
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
    
    def test_fails_cnsolidate(self):
        """ test that fails_cnsolidate() works correctly
        """
        
        # check when the variant is passed by CNSOLIDATE
        self.var.info["CNSOLIDATE"] = True
        self.assertFalse(self.var.fails_cnsolidate())
        
        # check when the variant hasn't been passed by CNSOLIDATE
        del self.var.info["CNSOLIDATE"]
        self.assertTrue(self.var.fails_cnsolidate())
    
    def test_fails_mad_ratio(self):
        """ test that fails_mad_ratio() works correctly
        """
        
        # check that var passes when MAD ratio > 15
        self.var.info["MEANLR2"] = "0.5"
        self.var.info["MADL2R"] = "0.02"
        self.assertFalse(self.var.fails_mad_ratio())
        
        # check that var passes when MAD ratio == 15
        self.var.info["MEANLR2"] = "0.3"
        self.var.info["MADL2R"] = "0.02"
        self.assertFalse(self.var.fails_mad_ratio())
        
        # check that var fails when MAD ratio < 15
        self.var.info["MEANLR2"] = "0.2"
        self.var.info["MADL2R"] = "0.02"
        self.assertTrue(self.var.fails_mad_ratio())
        
        # check that var fails when trying to divide by zero
        self.var.info["MEANLR2"] = "0.2"
        self.var.info["MADL2R"] = "0"
        self.assertTrue(self.var.fails_mad_ratio())
    
    def test_fails_wscore(self):
        """ test that fails_wscore() works correctly
        """
        
        # check that var passes when WSCORE > 0.4
        self.var.info["WSCORE"] = "0.5"
        self.assertFalse(self.var.fails_wscore())
        
        # check that var passes when WSCORE == 0.4
        self.var.info["WSCORE"] = "0.4"
        self.assertFalse(self.var.fails_wscore())
        
        # check that var fails when WSCORE > 0.4
        self.var.info["WSCORE"] = "0.399"
        self.assertTrue(self.var.fails_wscore())
    
    def test_fails_callp(self):
        """ test that fails_callp() works correctly
        """
        
        # check that var passes when CALLP < 0.01
        self.var.info["CALLP"] = "0.0"
        self.assertFalse(self.var.fails_callp())
        
        # check that var passes when CALLP == 0.01
        self.var.info["CALLP"] = "0.01"
        self.assertFalse(self.var.fails_callp())
        
        # check that var fails when CALLP > 0.01
        self.var.info["CALLP"] = "0.0101"
        self.assertTrue(self.var.fails_callp())
    
    def test_fails_commmon_forwards(self):
        """ test that fails_commmon_forwards() works correctly
        """
        
        # check that var passes when COMMONFORWARDS < 0.8
        self.var.info["COMMONFORWARDS"] = "0.0"
        self.assertFalse(self.var.fails_commmon_forwards())
        
        # check that var fails when COMMONFORWARDS == 0.8
        self.var.info["COMMONFORWARDS"] = "0.8"
        self.assertFalse(self.var.fails_commmon_forwards())
        
        # check that var fails when COMMONFORWARDS > 0.8
        self.var.info["COMMONFORWARDS"] = "0.9"
        self.assertTrue(self.var.fails_commmon_forwards())
    
    def test_fails_mean_lr2_dup(self):
        """ test that fails_mean_lr2() works correctly on duplications
        """
        
        # set the var as a duplication
        self.var.alt_allele = "<DUP>"
        self.var.set_genotype()
        
        # check that dup passes with MEANLR2 > 0.36
        self.var.info["MEANLR2"] = "0.4"
        self.assertFalse(self.var.fails_meanlr2())
        
        # check that dup passes with MEANLR2 == 0.36
        self.var.info["MEANLR2"] = "0.36"
        self.assertFalse(self.var.fails_meanlr2())
        
        # check that dup passes with MEANLR2 < 0.36
        self.var.info["MEANLR2"] = "0.359"
        self.assertTrue(self.var.fails_meanlr2())
    
    def test_fails_mean_lr2_del(self):
        """ test that fails_mean_lr2() works correctly on deletions
        """
        # set the var as a deletion
        self.var.alt_allele = "<DEL>"
        self.var.set_genotype()
        
        # check that del passes with MEANLR2 < -0.41
        self.var.info["MEANLR2"] = "-0.5"
        self.assertFalse(self.var.fails_meanlr2())
        
        # check that del passes with MEANLR2 == -0.41
        self.var.info["MEANLR2"] = "-0.41"
        self.assertFalse(self.var.fails_meanlr2())
        
        # check that del passes with MEANLR2 > -0.41
        self.var.info["MEANLR2"] = "-0.409"
        self.assertTrue(self.var.fails_meanlr2())
    
    def test_fails_no_exons(self):
        """ test that fails_no_exons() works correctly
        """
        
        self.var.info["NUMBEREXONS"] = "1"
        self.assertFalse(self.var.fails_no_exons())

        self.var.info["NUMBEREXONS"] = "0"
        self.assertTrue(self.var.fails_no_exons())


if __name__ == '__main__':
    unittest.main()


