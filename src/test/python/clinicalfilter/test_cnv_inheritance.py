""" unit testing of the Inheritance class
"""

import unittest
from clinicalfilter.ped import Family
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.inheritance import CNVInheritance
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
        self.known_genes = {"TEST": {"inh": {"Monoallelic": {"Loss of function"}}, "status": {"Confirmed DD Gene"}, "start": 5000, "end": 6000}}
        
        syndrome_regions = {("1", "1000", "2000"): 1}
        
        self.inh = CNVInheritance(self.trio, self.known_genes, "TEST", syndrome_regions)
    
    def create_cnv(self, gender, inh, cifer, chrom, pos):
        """ create a default variant
        """
        
        snp_id = "."
        ref = "A"
        alt = "<DUP>"
        filt = "PASS"
        
        # set up a SNV object, since SNV inherits VcfInfo
        var = CNV(chrom, pos, snp_id, ref, alt, filt)
        
        info = "HGNC=TEST;HGNC_ALL=TEST;END=16000000;SVLEN=5000"
        format_keys = "CIFER:INHERITANCE:DP"
        sample_values = cifer + ":" + inh + ":50"
        
        var.add_info(info)
        var.add_format(format_keys, sample_values)
        var.set_gender(gender)
        var.set_genotype()
        
        return var
    
    def create_variant(self, child_gender, chrom="1", position="15000000"):
        """ creates a TrioGenotypes variant
        """
        
        # generate a test variant
        child_var = self.create_cnv(child_gender, "unknown", "uncertain", chrom, position)
        mom_var = self.create_cnv("F", "unknown", "uncertain", chrom, position)
        dad_var = self.create_cnv("M", "unknown", "uncertain", chrom, position)
        
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
        
        cnv = self.create_variant("female")
        
        # check that paternally inherited CNVs that have affected fathers pass
        self.inh.trio.father.affected_status = "2"
        inh = ["paternal"]
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        
        # check for a CIFER derived annotation
        inh = ["paternal_inh"]
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        
        # check for a list of the VICAR and CIFER inheritance annotations
        inh = ["paternal", "paternal_inh"]
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        
        # check when one annotation says inherited, but the other doesn't.
        # Having one parentally inherited annotation is sufficient to classify
        # the CNV as inherited.
        inh = ["uknown", "paternal_inh"]
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        
        # check that paternally inherited CNVs without an affected father fail
        self.inh.trio.father.affected_status = "1"
        self.assertFalse(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        
        # check that maternally inherited CNVs without an affected mother fail
        inh = ["maternal"]
        self.inh.trio.father.affected_status = "1"
        self.assertFalse(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        
        # check that maternally inherited CNVs on chrX in male probands pass
        # regardless of the affected status of the mother.
        self.inh.trio.child.gender = "male"
        cnv.child.chrom = "X"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        
        # check that biparentally inherited CNVs pass if either parent is affected
        inh = ["biparental"]
        self.assertFalse(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        self.inh.trio.father.affected_status = "2"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        self.inh.trio.father.affected_status = "1"
        self.inh.trio.mother.affected_status = "2"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        
        # check that biparentally inherited CNVs pass if either parent is affected
        inh = ["inheritedDuo"]
        self.inh.trio.mother.affected_status = "1"
        self.assertFalse(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        self.inh.trio.father.affected_status = "2"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        self.inh.trio.father.affected_status = "1"
        self.inh.trio.mother.affected_status = "2"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
        
        # check that noninherited CNVs pass, regardless of parental affected status
        inh = ["deNovo"]
        self.inh.trio.mother.affected_status = "1"
        self.assertTrue(self.inh.inheritance_matches_parental_affected_status(cnv, inh))
    
    def test_passes_nonddg2p_filter(self):
        """ test that passes_nonddg2p_filter() works correctly
        """
        
        cnv = self.create_variant("female")
        cnv.child.genotype = "DUP"
        cnv.child.format["INHERITANCE"] = "deNovo"
        cnv.child.info["SVLEN"] = "1000001"
        
        # check that a sufficiently long de novo DUP passes
        self.assertTrue(self.inh.passes_nonddg2p_filter(cnv))
        
        # check that a insufficiently long de novo DUP fails
        cnv.child.info["SVLEN"] = "999999"
        self.assertFalse(self.inh.passes_nonddg2p_filter(cnv))
    
    def test_passes_gene_inheritance(self):
        """ test that passes_gene_inheritance() works correctly
        """
        
        gene = "TEST"
        inh = "Monoallelic"
        cnv = self.create_variant("female")
        
        cnv.chrom = "1"
        cnv.child.info["CNS"] = "3"
        cnv.child.gender = "M"
        cnv.child.genotype = "DUP"
        
        # check that a CNV with all the right characteristics passes
        # for mech in
        self.inh.known_genes[gene]["inh"][inh] = {"Increased gene dosage"}
        self.assertTrue(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        # check that a CNV in a gene with differing inheritance mechanism fails
        self.inh.known_genes[gene]["inh"][inh] = {"Loss of function"}
        self.assertTrue(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        # check that a DEL CNV requires a different mechanism
        cnv.child.genotype = "DEL"
        cnv.child.info["CNS"] = "0"
        self.assertTrue(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        self.inh.known_genes[gene]["inh"][inh] = {"Increased gene dosage"}
        self.assertTrue(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        # check that a CNV in a gene with "Uncertain" mechanism passes
        self.inh.known_genes[gene]["inh"][inh] = {"Uncertain"}
        self.assertTrue(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        # check that a CNV in a gene with "Uncertain" mechanism passes
        cnv.child.genotype = "DEL"
        cnv.child.info["CNS"] = "1"
        self.inh.known_genes[gene]["inh"][inh] = {"Uncertain"}
        self.assertTrue(self.inh.passes_gene_inheritance(cnv, gene, inh))
    
    def test_passes_gene_inheritance_biallelic(self):
        """ test that passes_gene_inheritance() works correctly for biallelic
        """
        
        gene = "TEST"
        inh = "Biallelic"
        cnv = self.create_variant("female")
        
        # check that a CNV with mismatched copy number fails
        cnv.child.genotype = "DEL"
        cnv.child.info["CNS"] = "3"
        self.inh.known_genes[gene]["inh"][inh] = {"Increased gene dosage"}
        self.assertFalse(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        # check that a CNV with correct copy number passes
        cnv.child.info["CNS"] = "0"
        self.inh.known_genes[gene]["inh"][inh] = {"Loss of function"}
        self.assertTrue(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        # check that a CNV with correct copy number, but wrong mechanism fails
        cnv.child.info["CNS"] = "0"
        self.inh.known_genes[gene]["inh"][inh] = {"Increased gene dosage"}
        self.assertFalse(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        # check that a CNV with mismatched copy number fails
        cnv.child.info["CNS"] = "1"
        self.inh.known_genes[gene]["inh"][inh] = {"Loss of function"}
        self.assertFalse(self.inh.passes_gene_inheritance(cnv, gene, inh))
    
    def test_passes_gene_inheritance_x_linked(self):
        """ test that passes_gene_inheritance() works correctly for X-linked dominant
        """
        
        gene = "TEST"
        inh = "X-linked dominant"
        cnv = self.create_variant("female")
        
        # check that a CNV with mismatched chrom fails
        cnv.child.genotype = "DUP"
        cnv.child.info["CNS"] = "3"
        self.inh.known_genes[gene]["inh"][inh] = {"Increased gene dosage"}
        self.assertFalse(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        # check that a CNV with correct chrom passes
        cnv.child.chrom = "X"
        cnv.child.info["CNS"] = "3"
        self.assertTrue(self.inh.passes_gene_inheritance(cnv, gene, inh))
    
    def test_passes_gene_inheritance_hemizygous(self):
        """ test that passes_gene_inheritance() works correctly for hemizygous
        """
        
        gene = "TEST"
        inh = "Hemizygous"
        cnv = self.create_variant("female")
        
        # check that female hemizygous CNV must be DUPs
        cnv.child.gender = "F"
        cnv.child.chrom = "X"
        cnv.child.info["CNS"] = "3"
        self.inh.known_genes[gene]["inh"][inh] = {"Increased gene dosage"}
        self.assertTrue(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        cnv.child.genotype = "DEL"
        cnv.child.info["CNS"] = "1"
        self.inh.known_genes[gene]["inh"][inh] = {"Loss of function"}
        self.assertFalse(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        # check that male hemizygous CNV can be either DEL or DUP
        cnv.child.gender = "M"
        self.assertTrue(self.inh.passes_gene_inheritance(cnv, gene, inh))
        
        cnv.child.genotype = "DUP"
        cnv.child.info["CNS"] = "3"
        self.inh.known_genes[gene]["inh"][inh] = {"Increased gene dosage"}
        self.assertTrue(self.inh.passes_gene_inheritance(cnv, gene, inh))
    
    def test_passes_gene_inheritance_unsupported_inh(self):
        """ test that passes_gene_inheritance() works correctly for unsupported inh
        """
        
        gene = "TEST"
        inh = "Mosaic"
        
        # check that non-supported inheritance modes fail, even if they
        # otherwise would
        cnv = self.create_variant("female")
        cnv.child.info["CNS"] = "1"
        self.inh.known_genes[gene]["inh"][inh] = {"Increased gene dosage"}
        self.assertFalse(self.inh.passes_gene_inheritance(cnv, gene, inh))
    
    def test_passes_ddg2p_filter(self):
        """ test if passes_ddg2p_filter() works correctly
        """
        
        gene_inh = {"TEST": {"inh": {"Monoallelic": \
            {"Increased gene dosage"}}, "status": \
            {"Confirmed DD Gene"}, "start": 5000, "end": 6000}}
        
        gene = "TEST"
        inh = "Monoallelic"
        cnv = self.create_variant("female")
        cnv.child.chrom = "1"
        cnv.child.info["CNS"] = "3"
        self.inh.known_genes[gene]["inh"][inh] = {"Increased gene dosage"}
        
        # check we don't get anything if there are no DDG2P genes
        self.inh.known_genes = {}
        self.assertFalse(self.inh.passes_ddg2p_filter(cnv))
        self.inh.known_genes = None
        self.assertFalse(self.inh.passes_ddg2p_filter(cnv))
        
        # check if the var passes when the inheritance mechanism, copy number
        # and chromosome are appropriate for the DDG2P gene
        self.inh.known_genes = gene_inh
        self.assertTrue(self.inh.passes_ddg2p_filter(cnv))
        
        # check that we only pass genes with exact matches from DDG2P
        self.inh.gene = "TEST1"
        self.assertFalse(self.inh.passes_ddg2p_filter(cnv))
        
        # check if the variant passes if the confirmed type is "Both DD and IF",
        # even if the variant wouldn't otherwise pass
        self.inh.gene = "TEST"
        self.inh.known_genes[gene]["status"] = {"Both DD and IF"}
        self.inh.known_genes[gene]["inh"][inh] = {"Loss of function"}
        self.assertTrue(self.inh.passes_ddg2p_filter(cnv))
        
        # fail on genes that don't have a robust confirmed status
        self.inh.known_genes[gene]["status"] = {"Possible DD Gene"}
        self.assertFalse(self.inh.passes_ddg2p_filter(cnv))
        
        # add another gene to the DDG2P dictionary
        self.inh.known_genes["TEST2"] = {"inh": {"Monoallelic": \
            {"Increased gene dosage"}}, "status": {"Both DD and IF"}}
        self.assertFalse(self.inh.passes_ddg2p_filter(cnv))
    
    def test_passes_ddg2p_filter_surrounding_disruptive_dup(self):
        pass
    
    def test_check_passes_intragenic_dup(self):
        """ test that passes_intragenic_dup() works correctly
        """
        
        gene = "TEST"
        inh = "Monoallelic"
        
        cnv = self.create_variant("female")
        
        # set parameters that will pass the function
        cnv.child.genotype = "DUP"
        cnv.position = "5200"
        cnv.child.info["END"] = "5800"
        
        gene_inh = {"TEST": {"inh": {"Monoallelic": \
            {"Loss of Function"}, "X-linked dominant": \
            {"Loss of Function"}}, "status": \
            {"Confirmed DD Gene"}, "start": 5000, "end": 6000}}
        
        self.inh.known_genes = gene_inh
            
        # check for values that pass the function
        self.assertTrue(self.inh.passes_intragenic_dup(cnv, gene, inh))
        
        # make a CNV that surrounds the entire gene, which will fail
        cnv.child.position = 4800
        cnv.child.info["END"] = "6200"
        self.assertFalse(self.inh.passes_intragenic_dup(cnv, gene, inh))
        
        # make a CNV where the gene exactly overlaps, which will fail
        cnv.child.position = 5000
        cnv.child.info["END"] = "6000"
        self.assertFalse(self.inh.passes_intragenic_dup(cnv, gene, inh))
        
        # make a CNV where the gene protrudes at the 5' end, which will pass
        cnv.child.position = 5200
        cnv.child.info["END"] = "6200"
        self.assertTrue(self.inh.passes_intragenic_dup(cnv, gene, inh))
        
        # make a CNV where the gene protrudes at the 3' end, which will pass
        cnv.child.position = 4800
        cnv.child.info["END"] = "5800"
        self.assertTrue(self.inh.passes_intragenic_dup(cnv, gene, inh))
        
        # check for correct response under different inheritance models
        self.assertTrue(self.inh.passes_intragenic_dup(cnv, gene, "X-linked dominant"))
        self.assertFalse(self.inh.passes_intragenic_dup(cnv, gene, "Biallelic"))
        
        # check that DEL CNVs fail
        cnv.child.genotype = "DEL"
        self.assertFalse(self.inh.passes_intragenic_dup(cnv, gene, inh))
    
    def test_check_cnv_region_overlap(self):
        """ test that check_cnv_region_overlap() works correctly
        """
        
        # set the variant key and copy number to known values
        cnv = self.create_variant("female")
        cnv.child.chrom = "1"
        cnv.child.position = 1000
        cnv.child.info["END"] = 2000
        cnv.child.info["CNS"] = "1"
        
        syndrome_regions = {("2", "5000", "6000"): 1, ("3", "8000", "9000"): 0}
        
        # check that if there aren't any overlapping regions, we return False
        self.assertFalse(self.inh.check_cnv_region_overlap(cnv, syndrome_regions))
        
        # check that when the region matches, but the chrom does not, we still
        # return False
        syndrome_regions[("2", "1000", "2000")] = "1"
        self.assertFalse(self.inh.check_cnv_region_overlap(cnv, syndrome_regions))
        
        # check that when the region and chrom overlap, but the copy number
        # does not, we still return False
        syndrome_regions[("1", "1000", "2000")] = "2"
        self.assertFalse(self.inh.check_cnv_region_overlap(cnv, syndrome_regions))
        
        # check that if the chrom, range and copy number overlap, and the
        # overlap region is sufficient, we return True
        syndrome_regions[("1", "1000", "2000")] = "1"
        self.assertTrue(self.inh.check_cnv_region_overlap(cnv, syndrome_regions))
    
    def test_has_enough_overlap(self):
        """ test that has_enough_overlap() works correctly
        """
        
        cnv_start = 1000
        cnv_end = 2000
        region_start = 1000
        region_end = 2000
        
        # check that CNV and syndrome region with 100% overlap, forwards and
        # backwards, pass
        self.assertTrue(self.inh.has_enough_overlap(cnv_start, cnv_end, region_start, region_end))
        
        # check that CNV with 50% overlap to the syndrome region passes
        cnv_end = 1500
        self.assertTrue(self.inh.has_enough_overlap(cnv_start, cnv_end, region_start, region_end))
        
        # check that CNV with < 50% overlap to the syndrome region fails
        cnv_end = 1499
        self.assertFalse(self.inh.has_enough_overlap(cnv_start, cnv_end, region_start, region_end))
        
        # check that syndrome region with 1% overlap to the CNV passes
        cnv_end = 2000
        region_end = 1010
        self.assertTrue(self.inh.has_enough_overlap(cnv_start, cnv_end, region_start, region_end))
        
        # check that 1 bp syndrome region which overlaps the CNV still passes
        region_end = 1001
        self.assertTrue(self.inh.has_enough_overlap(cnv_start, cnv_end, region_start, region_end))
    
    def test_check_compound_inheritance(self):
        """ test that check_compound_inheritance() works correctly
        """
        
        gene_inh = {"TEST": {"inh": {"Biallelic": \
            {"Increased gene dosage"}}, "status": \
            {"Confirmed DD Gene"}, "start": 5000, "end": 6000}}
        
        self.inh.known_genes = gene_inh
        cnv = self.create_variant("female")
        cnv.child.chrom = "1"
        cnv.child.info["CNS"] = "3"
        cnv.child.info["SVLEN"] = "1000001"
        
        # check that a standard CNV passes
        self.assertTrue(self.inh.check_compound_inheritance(cnv))
        
        # check that copy number = 0 passes
        cnv.child.info["CNS"] = "1"
        self.assertTrue(self.inh.check_compound_inheritance(cnv))
        
        # check that copy number = 0 fails
        cnv.child.info["CNS"] = "0"
        self.assertFalse(self.inh.check_compound_inheritance(cnv))
        
        # check that low SVLEN doesn't fail if the DDG2P route passes
        cnv.child.info["CNS"] = "1"
        cnv.child.info["SVLEN"] = "999999"
        self.assertTrue(self.inh.check_compound_inheritance(cnv))
        
        # check that low SVLEN combined with no DDG2P match fails
        del self.inh.known_genes["TEST"]["inh"]["Biallelic"]
        self.assertFalse(self.inh.check_compound_inheritance(cnv))
        
        # check that high SVLEN can overcome not having a DDG2P match
        cnv.child.info["SVLEN"] = "1000001"
        self.assertTrue(self.inh.check_compound_inheritance(cnv))
        
    def test_check_compound_inheritance_hemizygous(self):
        """ test that check_compound_inheritance() works correctly
        """
        
        gene_inh = {"TEST": {"inh": {"Hemizygous": \
            {"Increased gene dosage"}}, "status": \
            {"Confirmed DD Gene"}, "start": 5000, "end": 6000}}
        
        self.inh.known_genes = gene_inh
        cnv = self.create_variant("female")
        cnv.child.chrom = "X"
        cnv.child.info["CNS"] = "1"
        cnv.child.info["SVLEN"] = "500001"
        
        # check that a standard CNV passes
        self.assertTrue(self.inh.check_compound_inheritance(cnv))
        
        # check that low SVLEN doesn't fail if the DDG2P route passes
        cnv.child.info["SVLEN"] = "499999"
        self.assertTrue(self.inh.check_compound_inheritance(cnv))
        
        # check that low SVLEN combined with incorrect chrom fails
        cnv.child.chrom = "1"
        self.assertFalse(self.inh.check_compound_inheritance(cnv))
        
        # check that low SVLEN combined with incorrect sex fails
        cnv.child.chrom = "X"
        self.inh.trio.child.gender = "M"
        self.assertFalse(self.inh.check_compound_inheritance(cnv))
        
        # check that low SVLEN combined with incorrect copy number fails
        self.inh.trio.child.gender = "F"
        cnv.child.info["CNS"] = "3"
        self.assertFalse(self.inh.check_compound_inheritance(cnv))


if __name__ == '__main__':
    unittest.main()
