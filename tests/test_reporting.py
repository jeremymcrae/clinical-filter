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
import logging
import os
import datetime
import tempfile
import shutil
import gzip

import clinicalfilter
from clinicalfilter.ped import Family
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.variant.snv import SNV
from clinicalfilter.trio_genotypes import TrioGenotypes
from clinicalfilter.reporting import Report

logging.disable(logging.CRITICAL)

from tests.utils import create_snv, make_vcf_header

class TestReportPy(unittest.TestCase):
    """ test the Report class
    """
    
    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.mkdtemp()
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.temp_dir)
    
    def setUp(self):
        """ define a family and variant, and start the Allosomal class
        """
        
        # generate a test family
        child_gender = "F"
        mom_aff = "1"
        dad_aff = "1"
        
        self.trio = self.create_family(child_gender, mom_aff, dad_aff)
        
        # generate a test variant
        child = create_snv(child_gender, "0/1", chrom='X', pos=150,
            extra_info='MAX_AF=0.0005')
        mom = create_snv("F", "0/0", chrom='X', pos=150)
        dad = create_snv("M", "0/0", chrom='X', pos=150)
        
        self.variants = [TrioGenotypes('X', '150', child, mom, dad)]
        
        self.report = Report(None, None, None)
        self.report.family = self.trio
        SNV.set_populations(["AFR_AF", "AMR_AF", "ASN_AF", "DDD_AF",
            "EAS_AF", "ESP_AF", "EUR_AF", "MAX_AF", "SAS_AF", "UK10K_cohort_AF"])
    
    def create_family(self, child_gender, mom_aff, dad_aff):
        """ create a default family, with optional gender and parental statuses
        """
        
        fam = Family('test')
        fam.add_child('child', 'mother', 'father', child_gender, '2', 'child_vcf')
        fam.add_mother('mother', '0', '0', 'female', mom_aff, 'mother_vcf')
        fam.add_father('father', '0', '0', 'male', dad_aff, 'father_vcf')
        fam.set_child()
        
        return fam
    
    def test__save_tabular(self):
        ''' check that _save_tabular() works correctly
        '''
        
        temp = tempfile.NamedTemporaryFile(suffix='.txt', dir=self.temp_dir,
            delete=False)
        report = Report(temp.name, None, None)
        
        var = (self.variants[0], ["single_variant"], ["Monoallelic"], ["TEST"])
        var[0].child.format['GQ'] = 40
        report._save_tabular([var], self.trio)
        
        with open(temp.name, 'r') as handle:
            lines = handle.readlines()
        
        expected = ['proband\tsex\tchrom\tposition\tgene\t'
            'mutation_ID\ttranscript\tconsequence\tref/alt_alleles\tMAX_MAF\t'
            'inheritance\ttrio_genotype\tmom_aff\tdad_aff\tresult\tpp_dnm\t'
            'exac_allele_count\tGQ\thas_parents\tcnv_length\n',
            'child\tF\tX\t150\tTEST\tNA\tNA\t'
            'missense_variant\tA/G\t0.0005\tMonoallelic\t1/0/0\t1\t1\t'
            'single_variant\t0.99\tNA\t40\tTrue\tNA\n']
        
        self.assertEqual(lines, expected)
    
    def test__get_provenance(self):
        """ check that _get_provenance() works correctly
        """
        
        prov = ["checksum", "sample.calls.date.vcf.gz", "2014-01-01"]
        member = "proband"
        
        self.assertEqual(self.report._get_provenance(prov, member), \
            ["##UberVCF_proband_Id=sample\n", \
            "##UberVCF_proband_Checksum=checksum\n", \
            "##UberVCF_proband_Basename=sample.calls.date.vcf.gz\n", \
            "##UberVCF_proband_Date=2014-01-01\n"])
    
    def test__get_vcf_export_path(self):
        """ check that _get_vcf_export_path() works correctly
        """
        
        # use a folder to place the VCFG file in, which means we join the
        # proband ID to get a full path
        self.report.export_vcf = os.getcwd()
        self.assertEqual(self.report._get_vcf_export_path(),
            os.path.join(os.getcwd(), "child.vcf.gz"))
        
        # define an un-usable directory, to raise an error
        self.report.export_vcf = os.getcwd() + "asjhfgasjhfg"
        self.assertRaises(ValueError, self.report._get_vcf_export_path)
        
        # define a specific path for a VCF file, which is returned directly
        self.report.export_vcf = os.path.join(os.getcwd(), "sample_id.vcf.gz")
        self.assertEqual(self.report._get_vcf_export_path(), self.report.export_vcf)
    
    def test__make_vcf_header(self):
        """ check that _make_vcf_header() works correctly
        """
        
        # define the intial header lines
        header = make_vcf_header()
        
        # define the VCF provenances
        provenance = [("checksum", "proband.calls.date.vcf.gz", "2014-01-01"),
            ("checksum", "mother.calls.date.vcf.gz", "2014-01-02"),
            ("checksum", "father.calls.date.vcf.gz", "2014-01-03")]
        
        processed_header = ["##fileformat=VCFv4.1\n",
           '##fileDate=2014-01-01\n',
           "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
           '##INFO=<ID=ClinicalFilterType,Number=.,Type=String,'
                'Description="The type of clinical filter that passed this '
                'variant.">\n',
           '##INFO=<ID=ClinicalFilterGeneInheritance,Number=.,Type=String,'
                'Description="The inheritance mode (Monoallelic, Biallelic '
                'etc) under which the variant was found.">\n',
           '##INFO=<ID=ClinicalFilterReportableHGNC,Number=.,Type=String,'
                'Description="The HGNC symbol which the variant was identified '
                'as being reportable for.">\n',
           '##FORMAT=<ID=INHERITANCE_GENOTYPE,Number=.,Type=String,'
                'Description="The 012 coded genotypes for a trio (child, '
                'mother, father).">\n',
           '##FORMAT=<ID=INHERITANCE,Number=.,Type=String,Description="The '
                'inheritance of the variant in the trio (biparental, paternal, '
                'maternal, deNovo).">\n',
           "##ClinicalFilterRunDate={0}\n".format(datetime.date.today()),
           "##ClinicalFilterVersion={}\n".format(clinicalfilter.__version__),
           "##ClinicalFilterHistory=single_variant,compound_het\n",
           "##UberVCF_proband_Id=proband\n",
           "##UberVCF_proband_Checksum=checksum\n",
           "##UberVCF_proband_Basename=proband.calls.date.vcf.gz\n",
           "##UberVCF_proband_Date=2014-01-01\n",
           "##UberVCF_maternal_Id=mother\n",
           "##UberVCF_maternal_Checksum=checksum\n",
           "##UberVCF_maternal_Basename=mother.calls.date.vcf.gz\n",
           "##UberVCF_maternal_Date=2014-01-02\n",
           "##UberVCF_paternal_Id=father\n",
           "##UberVCF_paternal_Checksum=checksum\n",
           "##UberVCF_paternal_Basename=father.calls.date.vcf.gz\n",
           "##UberVCF_paternal_Date=2014-01-03\n",
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"]
        
        # check that the standard function returns the expected value. Note that
        # I haven't checked the output if self.known_genes_date is not None, nor
        # have I checked if the _clinicalFilterVersion is available
        self.assertEqual(self.report._make_vcf_header(header, provenance),
           processed_header)
    
    def test__get_parental_inheritance(self):
        """ check that _get_parental_inheritance() works correctly
        """
        
        var = self.variants[0]
        
        # check for the default genotypes
        self.assertEqual(self.report._get_parental_inheritance(var), "deNovo")
        
        # check when only the mother is non-ref
        var.mother.genotype = 1
        self.assertEqual(self.report._get_parental_inheritance(var), "maternal")
        
        # check when both parents are non-ref
        var.father.genotype = 1
        self.assertEqual(self.report._get_parental_inheritance(var), "biparental")
        
        # check when only the father is non-ref
        var.mother.genotype = 0
        self.assertEqual(self.report._get_parental_inheritance(var), "paternal")
        
        # check when the proband lacks parental information
        self.report.family.father = None
        self.report.family.mother = None
        self.assertEqual(self.report._get_parental_inheritance(var), "unknown")
    
    def test__get_vcf_lines(self):
        """ check that _get_vcf_lines() works correctly
        """
        
         # define the intial header lines
        header = make_vcf_header()
        
        # define the VCF provenances
        provenance = [("checksum", "proband.calls.date.vcf.gz", "2014-01-01"),
            ("checksum", "mother.calls.date.vcf.gz", "2014-01-02"),
            ("checksum", "father.calls.date.vcf.gz", "2014-01-03")]
        
        # define what the header will become
        vcf_lines = ["##fileformat=VCFv4.1\n",
           '##fileDate=2014-01-01\n',
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
           '##INFO=<ID=ClinicalFilterType,Number=.,Type=String,Description="The '
                'type of clinical filter that passed this variant.">\n',
           '##INFO=<ID=ClinicalFilterGeneInheritance,Number=.,Type=String,'
                'Description="The inheritance mode (Monoallelic, Biallelic etc) '
                'under which the variant was found.">\n',
           '##INFO=<ID=ClinicalFilterReportableHGNC,Number=.,Type=String,'
                'Description="The HGNC symbol which the variant was identified '
                'as being reportable for.">\n',
           '##FORMAT=<ID=INHERITANCE_GENOTYPE,Number=.,Type=String,'
                'Description="The 012 coded genotypes for a trio (child, '
                'mother, father).">\n',
           '##FORMAT=<ID=INHERITANCE,Number=.,Type=String,Description="'
                'The inheritance of the variant in the trio (biparental, '
                'paternal, maternal, deNovo).">\n',
           "##ClinicalFilterRunDate={0}\n".format(datetime.date.today()),
           "##ClinicalFilterVersion={}\n".format(clinicalfilter.__version__),
           "##ClinicalFilterHistory=single_variant,compound_het\n",
           "##UberVCF_proband_Id=proband\n",
           "##UberVCF_proband_Checksum=checksum\n",
           "##UberVCF_proband_Basename=proband.calls.date.vcf.gz\n",
           "##UberVCF_proband_Date=2014-01-01\n",
           "##UberVCF_maternal_Id=mother\n",
           "##UberVCF_maternal_Checksum=checksum\n",
           "##UberVCF_maternal_Basename=mother.calls.date.vcf.gz\n",
           "##UberVCF_maternal_Date=2014-01-02\n",
           "##UberVCF_paternal_Id=father\n",
           "##UberVCF_paternal_Checksum=checksum\n",
           "##UberVCF_paternal_Basename=father.calls.date.vcf.gz\n",
           "##UberVCF_paternal_Date=2014-01-03\n",
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"]
        
        # define what the default variant vcf line will become
        line = ['X\t150\t.\tA\tG\t50\tPASS\tCQ=missense_variant;'
            'ClinicalFilterGeneInheritance=Monoallelic;'
            'ClinicalFilterReportableHGNC=TEST;ClinicalFilterType=single_variant;'
            'DENOVO-SNP;HGNC=TEST;MAX_AF=0.0005\tGT:DP:INHERITANCE:'
            'INHERITANCE_GENOTYPE\t0/1:50:deNovo:1,0,0\n']
        
        # check that a list of one variant produces the correct VCF output. Note
        # that we haven't checked against CNVs, which can change the
        # INHERITANCE_GENOTYPE flag, nor have we tested a larger list of variants
        var = (self.variants[0], ["single_variant"], ["Monoallelic"], ["TEST"])
        var[0].child.add_vcf_line(['X', '150', '.', 'A', 'G', '50',
            'PASS', 'HGNC=TEST;CQ=missense_variant;EUR_AF=0.0005',
            'GT:DP', '0/1:50'])
        
        self.assertEqual(self.report._get_vcf_lines([var], header, provenance), vcf_lines + line)
    
    def test__get_output_line(self):
        """ check that _get_output_line() works correctly
        """
        
        var = (self.variants[0], ["single_variant"], ["Monoallelic"], ["TEST"])
        
        # check the output for the default variant
        expected = "child\tF\tX\t150\tTEST\tNA\tNA\tmissense_variant\t" \
            "A/G\t0.0005\tMonoallelic\t1/0/0\t1\t1\tsingle_variant\t0.99\tNA\tNA\tTrue\tNA\n"
        self.assertEqual(self.report._get_output_line(var, self.trio), expected)
        
        # introduce additional info for the output line parsing, check the line
        # that is returned is expected
        var[0].child.info["PolyPhen"] = "probably_damaging(0.99)"
        var[0].child.info["SIFT"] = "deleterious(0)"
        var[0].child.info["ENST"] = "ENST00X"
        expected = "child\tF\tX\t150\tTEST\tNA\tENST00X\t" \
            "missense_variant,PolyPhen=probably_damaging(0.99)," \
            "SIFT=deleterious(0)\tA/G\t0.0005\tMonoallelic\t1/0/0\t1\t1\t" \
            "single_variant\t0.99\tNA\tNA\tTrue\tNA\n"
        self.assertEqual(self.report._get_output_line(var, self.trio), expected)
    
    def test__write_vcf(self):
        ''' check that _write_vcf() works correctly
        '''
        
        path = tempfile.NamedTemporaryFile(suffix='.vcf.gz', dir=self.temp_dir,
            delete=False)
        lines = make_vcf_header() +  ['X\t150\t.\tA\tG\t50\tPASS\tHGNC=TEST;'
            'CQ=missense_variant;random_tag;EUR_AF=0.0005;'
            'ClinicalFilterGeneInheritance=Monoallelic;'
            'ClinicalFilterType=single_variant;'
            'ClinicalFilterReportableHGNC=TEST\tGT:DP:INHERITANCE:'
            'INHERITANCE_GENOTYPE\t0/1:50:deNovo:1,0,0\n']
        
        self.report._write_vcf(path.name, lines)
        
        with gzip.open(path.name, 'r') as handle:
            vcf = [ x.decode() for x in handle ]
            self.assertEqual(lines, vcf)
        
