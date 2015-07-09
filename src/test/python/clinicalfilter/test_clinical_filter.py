""" unit testing of the ClinicalFilter class
"""

import unittest

from clinical_filter import ClinicalFilter
from clinicalfilter.ped import Family
from clinicalfilter.load_options import LoadOptions
from clinicalfilter.variant.snv import SNV
from clinicalfilter.trio_genotypes import TrioGenotypes


class TestClinicalFilterPy(unittest.TestCase):
    """ test the ClinicalFilter class
    """
    
    def setUp(self):
        """ create a default ClinicalFilter object to test
        """
        
        class opts:
            ped = None
            child, mother, father = None, None, None
            gender = None
            mom_aff, dad_aff = None, None
            regions = None
            genes, genes_date = None, None
            alternate_ids = None
            output, export_vcf = None, None
            loglevel = "debug"
            debug_chrom, debug_pos = None, None
            pp_filter = 0.9
        
        self.finder = ClinicalFilter(opts)
    
    def create_snv(self, sex, genotype, cq="missense_variant", hgnc="TEST", chrom="1"):
        """ create a default variant
        """
        
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "G"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = SNV(chrom, pos, snp_id, ref, alt, filt)
        
        info = "HGNC={0};CQ={1};DENOVO-SNP;PP_DNM=0.99".format(hgnc, cq)
        keys = "GT:DP:TEAM29_FILTER:PP_DNM"
        values = genotype + ":50:PASS:0.99"
        
        var.add_info(info)
        var.add_format(keys, values)
        var.set_gender(sex)
        var.set_genotype()
        
        return var
    
    def create_trio_variant(self, child_gender, cq, hgnc, chrom="1"):
        """ create a default TrioGenotypes variant
        """
        
        # generate a test variant
        child_var = self.create_snv(child_gender, "0/1", cq, hgnc, chrom)
        mom_var = self.create_snv("F", "0/0", cq, hgnc, chrom)
        dad_var = self.create_snv("M", "0/0", cq, hgnc, chrom)
        
        var = TrioGenotypes(child_var)
        var.add_mother_variant(mom_var)
        var.add_father_variant(dad_var)
        
        return var
    
    def test_create_gene_dict(self):
        """ test that create_gene_dict works correctly
        """
        
        # create variants that share genes, or not
        snv1 = self.create_trio_variant("F", "missense_variant|missense_variant", "TEST1|TEST2")
        snv2 = self.create_trio_variant("F", "missense_variant", "TEST1")
        snv3 = self.create_trio_variant("F", "missense_variant", "OTHER1")
        
        # the variants that share a gene should be grouped in lists indexed by
        # the gene key
        self.assertEqual(self.finder.create_gene_dict([snv1, snv2, snv3]),
            {"TEST1": [snv1, snv2], "TEST2": [snv1], "OTHER1": [snv3]})
    
    def test_find_variants(self):
        """ test that find_variants() works correctly
        """
        
        # define the trio, so that we can know whether the parents are affected.
        # The child also needs to be included and set, so that we can get the
        # child ID for logging purposes.
        self.finder.family = Family("famID")
        self.finder.family.add_child("child_id", "/vcf/path", "2", "F")
        self.finder.family.add_father("dad_id", "/vcf/path", "1", "M")
        self.finder.family.add_mother("mom_id", "/vcf/path", "1", "F")
        self.finder.family.set_child()
        
        # create variants that cover various scenarios
        snv1 = self.create_trio_variant("F", "missense_variant|missense_variant", "TEST1|TEST2")
        snv2 = self.create_trio_variant("F", "missense_variant|synonymous_variant", "OTHER1|OTHER2")
        snv3 = self.create_trio_variant("F", "missense_variant", "")
        snv4 = self.create_trio_variant("F", "missense_variant", "TESTX", chrom="X")
        
        self.finder.known_genes = {"TEST1": {"inh": ["Monoallelic"]},
            "OTHER1": {"inh": ["Monoallelic"]},
            "OTHER2": {"inh": ["Monoallelic"]},
            "TESTX": {"inh": ["X-linked dominant"]}}
        
        # check the simplest case, a variant in a known gene
        self.assertEqual(self.finder.find_variants([snv1], "TEST1"),
            [(snv1, "single_variant", "Monoallelic", ["TEST1"])])
        
        # check that a gene not in a known gene does not pass
        self.assertEqual(self.finder.find_variants([snv1], "TEST2"), [])
        
        # check a variant where the gene is known, but the consequence for that
        # gene is not functional, does not pass
        self.assertEqual(self.finder.find_variants([snv2], "OTHER2"), [])
        
        # check that intergenic variants (which lack HGNC symbols) do not pass
        self.assertEqual(self.finder.find_variants([snv3], None), [])
        
        # check that a variant on chrX passes through the allosomal instance
        self.assertEqual(self.finder.find_variants([snv4], "TESTX"),
            [(snv4, "single_variant", "X-linked dominant", ["TESTX"])])
        
        # remove the known genes, so that the variants in unknown genes pass
        self.finder.known_genes = None
        self.assertEqual(self.finder.find_variants([snv1], "TEST2"),
            [(snv1, "single_variant", "Monoallelic", ["TEST2"])])
        
        # but variants without gene symbols still are excluded
        self.assertEqual(self.finder.find_variants([snv3], None), [])
    
    def test_exclude_duplicates(self):
        """ test that exclude duplicates works correctly
        """
        
        # create a variant that is within two genes
        snv1 = self.create_trio_variant("F", "missense_variant|missense_variant", "TEST1|TEST2")
        
        # two variants that lie in different genes on different chromosomes
        # should not be merged
        snv2 = self.create_trio_variant("F", "missense_variant", "OTHER1", chrom="2")
        variants = [(snv1, "single_variant", "Monoallelic", ["TEST1"]),
            ((snv2, "single_variant", "Monoallelic", ["OTHER1"]))]
        self.assertEqual(sorted(self.finder.exclude_duplicates(variants)), sorted(variants))
        
        # create a list of variant tuples that passed filtering for two
        # different gene symbols
        variants = [(snv1, "single_variant", "Monoallelic", ["TEST1"]),
            ((snv1, "compound_het", "Biallelic", ["TEST1"])),
            ((snv1, "compound_het", "Biallelic", ["TEST1"]))]
        self.assertEqual(self.finder.exclude_duplicates(variants),
            [(snv1, "single_variant,compound_het", "Monoallelic,Biallelic", ["TEST1"])])
        
        # create a list of variant tuples that passed filtering for two
        # different gene symbols
        variants = [(snv1, "single_variant", "Monoallelic", ["TEST1"]),
            ((snv1, "single_variant", "Monoallelic", ["TEST2"]))]
        
        # the same variant passing for two gene symbols should be collapsed
        # into a single entry, where the entry contains a list ofall the gene
        # symbols
        self.assertEqual(self.finder.exclude_duplicates(variants),
            [(snv1, "single_variant", "Monoallelic", ["TEST1", "TEST2"])])
        
        
        
        
