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
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO
import sys

from clinicalfilter.variant.snv import SNV

class TestVariantSnvPy(unittest.TestCase):
    """ unit testing of the SNV class
    """
    
    pops = ["AFR_AF", "AMR_AF", "ASN_AF", "DDD_AF", "EAS_AF", "ESP_AF",
        "EUR_AF", "MAX_AF", "SAS_AF", "UK10K_cohort_AF"]
    SNV.set_populations(pops)
    
    def setUp(self):
        """ define a default VcfInfo object
        """
        
        chrom = "1"
        pos = "15000000"
        snp_id = "."
        ref = "A"
        alt = "G"
        filt = "PASS"
        
        info = "HGNC=ATRX;CQ=missense_variant;random_tag"
        self.keys = "GT:DP"
        self.values = "0/1:50"
        
        # set up a SNV object, since SNV inherits VcfInfo
        self.var = SNV(chrom, pos, snp_id, ref, alt, filt, info=info,
            format=self.keys, sample=self.values)
    
    def test_get_key(self):
        """ tests that get_key() operates correctly
        """
        
        # make sure the chrom and position are correct
        self.var.chrom = "1"
        self.var.position = "15000000"
        
        self.assertEqual(self.var.get_key(), ("1", "15000000"))
        
        # and make sure the chrom and position are correct if we change them
        self.var.chrom = "22"
        self.var.position = "123456789"
        self.assertEqual(self.var.get_key(), ("22", "123456789"))
    
    def test_convert_genotype(self):
        """ test that genotypes convert from two char to single char
        """
        
        genotypes = [("0/0", 0), ("0/1", 1), ("1/0", 1), ("1/1", 2), \
            ("1/2", 1), ("2/1", 1), ("0/2", 1), ("2/0", 1), ("2/2", 2)]
        
        # run thorugh all the legit genotype codes
        for geno in genotypes:
            genotype = geno[0]
            result = geno[1]
            self.assertEqual(self.var.convert_genotype(genotype), result)
         
        # Raise error when converting single character genotype
        with self.assertRaises(ValueError):
            self.var.convert_genotype("0")
          
        # raise error when converting unknown genotype
        with self.assertRaises(AssertionError):
            self.var.convert_genotype("a/a")
            
        # also include other genotype format posibilities. None of these are
        # used, but since they aren't explicitly forbidden, make sure they work
        
        # check two character strings
        self.assertEqual(self.var.convert_genotype("12|34"), 1)
        self.assertEqual(self.var.convert_genotype("99|99"), 2)
    
    def test_set_genotype_autosomal(self):
        """ test that set_genotype() operates correctly
        """
        
        self.var.add_format(self.keys, self.values)
        self.var._set_gender("male")
        
        genotypes = [("0/0", 0), ("0/1", 1), ("1/1", 2)]
        
        for geno in genotypes:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.get_genotype(), result)
        
        # remove the format attribute, so we can raise an error
        self.var.format = None
        with self.assertRaises(ValueError):
            self.var.set_genotype()
    
    def test_set_genotype_allosomal_male(self):
        """ test that set_genotype() operates correctly for the male X chrom
        """
        
        self.var.add_format(self.keys, self.values)
        self.var.chrom = "X"
        self.var._set_gender("male")
        
        genotypes = [("0/0", 0), ("1/1", 2)]
        
        for geno in genotypes:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.get_genotype(), result)
        
        # check that we raise an error for X chrom hets
        genotypes = ["0/1", "1/0"]
        for genotype in genotypes:
            self.var.format["GT"] = genotype
            with self.assertRaises(ValueError):
                self.var.set_genotype()
    
    def test_set_genotype_allosomal_female(self):
        """ test that set_genotype() operates correctly for the female X chrom
        """
        
        self.var.add_format(self.keys, self.values)
        self.var.chrom = "X"
        self.var._set_gender("female")
        
        genotypes = [("0/0", 0), ("0/1", 1), ("1/1", 2)]
        
        for geno in genotypes:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.get_genotype(), result)
    
    def test_is_het_autosomal(self):
        """ tests that is_het() operates correctly for automsal chromosomes
        """
        
        self.var.add_format(self.keys, self.values)
        self.var._set_gender("male")
        
        het = [("0/0", False), ("0/1", True), ("1/1", False)]
        
        for geno in het:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.is_het(), result)
     
    def test_is_hom_alt_autosomal(self):
        """ tests that is_hom_alt() operates correctly for automsal chromosomes
        """
        
        self.var.add_format(self.keys, self.values)
        self.var._set_gender("male")
        
        hom_alt = [("0/0", False), ("0/1", False), ("1/1", True)]
        
        for geno in hom_alt:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.is_hom_alt(), result)
            
    def test_is_hom_ref_autosomal(self):
        """ tests that is_hom_ref() operates correctly for automsal chromosomes
        """
        
        self.var.add_format(self.keys, self.values)
        self.var._set_gender("male")
        
        hom_ref = [("0/0", True), ("0/1", False), ("1/1", False)]
        
        for geno in hom_ref:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.is_hom_ref(), result)
        
    def test_is_not_ref_autosomal(self):
        """ tests that is_not_ref() operates correctly for automsal chromosomes
        """
        
        self.var.add_format(self.keys, self.values)
        self.var._set_gender("male")
        
        not_ref = [("0/0", False), ("0/1", True), ("1/1", True)]
        
        for geno in not_ref:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.is_not_ref(), result)
    
    def test_is_not_alt_autosomal(self):
        """ tests that is_not_ref() operates correctly for automsal chromosomes
        """
        
        self.var.add_format(self.keys, self.values)
        self.var._set_gender("male")
        
        not_alt = [("0/0", True), ("0/1", True), ("1/1", False)]
        
        for geno in not_alt:
            genotype = geno[0]
            result = geno[1]
            
            self.var.format["GT"] = genotype
            self.var.set_genotype()
            self.assertEqual(self.var.is_not_alt(), result)
    
    def test_passes_default_filters(self):
        """ test that different variants pass or fail the VcfInfo filters
        """
        
        # check that a default variant passes the filters
        self.assertTrue(self.var.passes_filters())
    
    def test_passes_alternate_filter_string(self):
        """ test that the alternate permitted FILTER string also passes
        """
        
        # check that the alternate FILTER value passes
        self.var.filter = "."
        self.assertTrue(self.var.passes_filters())
        
        self.var.filter = "FAIL"
        self.assertFalse(self.var.passes_filters())
        
        # check that low VQSLOD on its own will pass the variant
        self.var.filter = "LOW_VQSLOD"
        self.assertTrue(self.var.passes_filters())
        
        # check that low VQSLOD in a de novo will still pass
        self.var.filter = "LOW_VQSLOD"
        self.var.info["DENOVO-SNP"] = True
        self.assertTrue(self.var.passes_filters())
    
    def test_passes_filters_low_maf(self):
        """ tests that low MAF values pass the filters
        """
        
        # check that low MAF values pass the filters
        for pop in self.pops:
            self.var.info[pop] = "0.001"
            self.assertTrue(self.var.passes_filters())
            
            # and check that MAF on the threshold still pass
            self.var.info[pop] = "0.005"
            self.assertTrue(self.var.passes_filters())
    
    def test_out_of_range_maf(self):
        """ check that MAF outside 0-1 still pass or fail correctly
        """
        
        self.var.info["AFR_AF"] = "-1"
        self.assertTrue(self.var.passes_filters())
      
        self.var.info["AFR_AF"] = "100"
        self.assertFalse(self.var.passes_filters())
      
    def test_fails_filters_high_maf(self):
        """ test that variants with high MAF fail the filtering
        """
        
        # check th
        for pop in self.pops:
            var = self.var
            var.info[pop] = "0.0101"
            self.assertFalse(var.passes_filters())
    
    def test_passes_consequence_filter(self):
        """ check all the consequence values that should pass
        """
        
        vep_passing = ["transcript_ablation", "splice_donor_variant", \
            "splice_acceptor_variant", "frameshift_variant", \
            "initiator_codon_variant", "inframe_insertion", "inframe_deletion",\
            "missense_variant", "transcript_amplification", "stop_gained",\
            "stop_lost"]
        
        # check all the passing consequences
        for cq in vep_passing:
            self.var.consequence = [cq]
            self.assertTrue(self.var.passes_filters())
            
    def test_fails_consequence_filter(self):
        """ check all the consequence values that should fail
        """
        
        vep_failing = ["splice_region_variant", \
            "incomplete_terminal_codon_variant", "synonymous_variant", \
            "stop_retained_variant", "mature_miRNA_variant", \
            "5_prime_UTR_variant", "3_prime_UTR_variant", \
            "non_coding_exon_variant", "nc_transcript_variant", \
            "intron_variant", "NMD_transcript_variant", \
            "upstream_gene_variant", "downstream_gene_variant", \
            "TFBS_ablation", "TFBS_amplification", "TF_binding_site_variant", \
            "regulatory_region_variant", "regulatory_region_ablation", \
            "regulatory_region_amplification", "feature_elongation", \
            "feature_truncation", "intergenic_variant", "coding_sequence_variant"]
        
        # check all the failing consequences
        for cq in vep_failing:
            self.var.consequence = [cq]
            self.assertFalse(self.var.passes_filters())
    
    def test_passes_filters_with_debug(self):
        """ check that passes_filters_with_debug() generates a failure message
        """
        
        # make a variant that will fail the filtering, and set the site for
        # debugging
        self.var.info["AFR_AF"] = "0.05"
        self.var.debug_pos = self.var.get_position()
        
        # get ready to capture the output from a print function
        out = StringIO()
        sys.stdout = out
        
        # check that the variant fails (and secondarily prints the failure mode)
        self.assertFalse(self.var.passes_filters_with_debug())
        output = out.getvalue().strip()
        
        # check that the message about why the variant failed filtering is correct
        self.assertEqual(output, "failed MAF: 0.05")
        
        # reset the standard out, so that we can observe other print statements
        sys.stdout = sys.__stdout__
    
