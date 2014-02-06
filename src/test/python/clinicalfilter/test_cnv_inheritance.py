""" unit testing of the Inheritance class
"""

import unittest
from clinicalfilter.ped import Family
from clinicalfilter.ped import Person
from clinicalfilter.variant import Variant
from clinicalfilter.variant_cnv import CNV
from clinicalfilter.inheritance import CNVInheritance
from clinicalfilter.vcf_info import VcfInfo
from clinicalfilter.trio_genotypes import TrioGenotypes


class TestCNVInheritancePy(unittest.TestCase):
    """ test the Inheritance class
    """
    
    def setUp(self):
        """ define a family and variant, and start the Inheritance class
        """
        
        # generate a test family
        child_gender = "F"
        mom_aff = "1"
        dad_aff = "1"
        
        self.trio = self.create_family(child_gender, mom_aff, dad_aff)
        
        # generate list of variants
        self.variant = self.create_variant(child_gender)
        
        # make sure we've got known genes data
        self.known_genes = {"TEST": {"inheritance": {"Monoallelic": {"Loss of function"}}, "confirmed_status": {"Confirmed DD Gene"}}}
        
        self.inh = CNVInheritance(self.variant, self.trio, self.known_genes)
    
    def create_cnv(self, gender, inh, chrom, pos):
        """ create a default variant
        """
        
        snp_id = "."
        ref = "A"
        alt = "<DUP>"
        qual = "50"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = CNV(chrom, pos, snp_id, ref, alt, qual, filt)
        
        tags = {"gene": ["HGNC", "VGN", "GN"], "consequence": ["VCQ", "CQ"]}
        
        info = "HGNC=TEST;HGNC_ALL=TEST;END=16000000;SVLEN=5000"
        format_keys = "inheritance:DP"
        sample_values = inh + ":50"
        
        var.add_info(info, tags)
        var.add_format(format_keys, sample_values)
        var.set_gender(gender)
        var.set_genotype()
        
        return var
    
    def create_variant(self, child_gender, chrom="1", position="15000000"):
        """ creates a TrioGenotypes variant
        """
        
        # generate a test variant
        child_var = self.create_cnv(child_gender, "0/1", chrom, position)
        mom_var = self.create_cnv("F", "0/0", chrom, position)
        dad_var = self.create_cnv("M", "0/0", chrom, position)
        
        var = TrioGenotypes(child_var)
        var.add_mother_variant(mom_var)
        var.add_father_variant(dad_var)
        
        return var
    
    def create_family(self, child_gender, mom_aff, dad_aff):
        """ create a default family, with optional gender and parental statuses
        """
        
        fam = Family("test")
        fam.add_child("child", "child_vcf", "2", child_gender)
        fam.add_mother("mother", "mother_vcf", mom_aff, "2")
        fam.add_father("father", "father_vcf", dad_aff, "1")
        fam.set_child()
        
        return fam
    
    def test_inheritance_matches_parental_affected_status(self):
        """ test that inheritance_matches_parental_affected_status() works correctly
        """
        
        # check that paternally inherited CNVs that have affected fathers pass
        self.inh.trio.father.affected_status = "2"
        inh = "paternal"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(inh))
        
        # check that paternally inherited CNVs without an affected father fail
        self.inh.trio.father.affected_status = "1"
        self.assertFalse(self.inh.inheritance_matches_parental_affected_status(inh))
        
        # check that maternally inherited CNVs without an affected mother fail
        inh = "maternal"
        self.inh.trio.father.affected_status = "1"
        self.assertFalse(self.inh.inheritance_matches_parental_affected_status(inh))
        
        # check that biparentally inherited CNVs pass if either parent is affected
        inh = "biparental"
        self.assertFalse(self.inh.inheritance_matches_parental_affected_status(inh))
        self.inh.trio.father.affected_status = "2"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(inh))
        self.inh.trio.father.affected_status = "1"
        self.inh.trio.mother.affected_status = "2"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(inh))
        
        # check that biparentally inherited CNVs pass if either parent is affected
        inh = "inheritedDuo"
        self.inh.trio.mother.affected_status = "1"
        self.assertFalse(self.inh.inheritance_matches_parental_affected_status(inh))
        self.inh.trio.father.affected_status = "2"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(inh))
        self.inh.trio.father.affected_status = "1"
        self.inh.trio.mother.affected_status = "2"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(inh))
        
        # check that noninherited CNVs pass, regardless of parental affected status
        inh = "deNovo"
        self.inh.trio.mother.affected_status = "1"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(inh))
    
    def test_passes_nonddg2p_filter(self):
        """ test that passes_nonddg2p_filter() works correctly
        """
        
        self.inh.variant.child.genotype = "DUP"
        self.inh.variant.child.format["INHERITANCE"] = "deNovo"
        self.inh.variant.child.info["SVLEN"] = "250001"
        
        # check that a sufficiently long de novo DUP passes
        self.assertTrue(self.inh.passes_nonddg2p_filter())
        
        # check that a insufficiently long de novo DUP fails
        self.inh.variant.child.info["SVLEN"] = "249999"
        self.assertFalse(self.inh.passes_nonddg2p_filter())
        
        # check that a sufficiently long de novo DEL passes
        self.inh.variant.child.genotype = "DEL"
        self.assertTrue(self.inh.passes_nonddg2p_filter())
        
        # check that a insufficiently long de novo DEL fails
        self.inh.variant.child.info["SVLEN"] = "99999"
        self.assertFalse(self.inh.passes_nonddg2p_filter())
        
        # check that paternally inherited CNVs that have affected fathers pass
        self.inh.variant.child.info["SVLEN"] = "600000"
        self.inh.trio.father.affected_status = "2"
        self.inh.variant.child.format["INHERITANCE"] = "paternal"
        self.assertTrue(self.inh.passes_nonddg2p_filter())
        
        # check that non-inherited CNVs pass a more stringent length
        self.inh.variant.child.format["INHERITANCE"] = "inconclusive"
        self.assertTrue(self.inh.passes_nonddg2p_filter())
        self.inh.variant.child.info["SVLEN"] = "499999"
        self.assertFalse(self.inh.passes_nonddg2p_filter())
    
    def test_passes_gene_inheritance(self):
        """ test that passes_gene_inheritance() works correctly
        """
        
        gene = "TEST"
        inh = "Monoallelic"
        
        self.inh.variant.chrom = "1"
        self.inh.variant.child.genotype = "DUP"
        self.inh.variant.child.info["CNS"] = "3"
        self.inh.known_genes[gene]["inheritance"][inh] = {"Increased gene dosage"}
        
        self.assertTrue(self.inh.passes_gene_inheritance(gene, inh))
        
        
        
        
        
        
        
        
        


if __name__ == '__main__':
    unittest.main()

