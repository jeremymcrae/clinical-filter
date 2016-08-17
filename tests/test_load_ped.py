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
import tempfile

from clinicalfilter.ped import load_ped_file, load_families, Person, Family

class TestLoadPed(unittest.TestCase):
    """ unit testing of loading ped files
    """
    
    def setUp(self):
        """ define and write a temporary ped file
        """
        
        self.tempfile = tempfile.NamedTemporaryFile(mode="w")
        self.path = self.tempfile.name
        
        self.tempfile.write("fam_ID   proband   dad   mom   F  2  /path/to/proband_vcf.gz\n")
        self.tempfile.write("fam_ID   dad       0     0     M  1  /path/to/dad_vcf.gz\n")
        self.tempfile.write("fam_ID   mom       0     0     F  1  /path/to/mom_vcf.gz\n")
        self.tempfile.flush()
    
    def test_load_ped_file_single_family(self):
        """ check that we correctly parse a ped file with a single trio
        """
        
        # load all the components from the file
        mothers, fathers, children, affected, sex, vcfs = load_ped_file(self.path)
        
        # check that they all match
        self.assertEqual(mothers, {'mom': '0', 'proband': 'mom', 'dad': '0'})
        self.assertEqual(fathers, {'mom': '0', 'proband': 'dad', 'dad': '0'})
        self.assertEqual(children, {'proband': 'fam_ID'})
        self.assertEqual(affected, {'mom': '1', 'proband': '2', 'dad': '1'})
        self.assertEqual(sex, {'mom': 'F', 'proband': 'F', 'dad': 'M'})
        self.assertEqual(vcfs, {'mom': '/path/to/mom_vcf.gz', \
            'proband': '/path/to/proband_vcf.gz', \
            'dad': '/path/to/dad_vcf.gz'})
    
    def test_load_ped_file_multiple_sibs(self):
        """ check that we correctly parse a ped file with multiple siblings
        """
        
        # add an extra sibling
        self.tempfile.write("fam_ID  sib   dad  mom  F  2  /path/to/sib_vcf.gz\n")
        self.tempfile.flush()
        
        # load all the components from the file
        mothers, fathers, children, affected, sex, vcfs = load_ped_file(self.path)
        
        # check that they all match
        self.assertEqual(mothers, {'mom': '0', 'proband': 'mom', 'dad': '0', 'sib': 'mom'})
        self.assertEqual(fathers, {'mom': '0', 'proband': 'dad', 'dad': '0', 'sib': 'dad'})
        self.assertEqual(children, {'proband': 'fam_ID', "sib": "fam_ID"})
        self.assertEqual(affected, {'mom': '1', 'proband': '2', 'dad': '1', 'sib': '2'})
        self.assertEqual(sex, {'mom': 'F', 'proband': 'F', 'dad': 'M', 'sib': 'F'})
        self.assertEqual(vcfs, {'mom': '/path/to/mom_vcf.gz',
            'proband': '/path/to/proband_vcf.gz',
            'dad': '/path/to/dad_vcf.gz',
            'sib': '/path/to/sib_vcf.gz'})
    
    def test_load_ped_file_multiple_families(self):
        """ check that we correctly parse a ped file with multiple families
        """
        
        # add an extra family, with multiple sibs
        self.tempfile.write("fam_ID2  proband2  dad2  mom2  F  2  /path/to/proband2_vcf.gz\n")
        self.tempfile.write("fam_ID2  dad2      0     0     M  1  /path/to/dad2_vcf.gz\n")
        self.tempfile.write("fam_ID2  mom2      0     0     F  1  /path/to/mom2_vcf.gz\n")
        self.tempfile.write("fam_ID2  sib       dad2  mom2  F  2  /path/to/sib_vcf.gz\n")
        self.tempfile.flush()
        
        # load all the components from the file
        mothers, fathers, children, affected, sex, vcfs = load_ped_file(self.path)
        
        # check that they all match
        self.assertEqual(mothers, {'mom': '0', 'proband': 'mom', 'dad': '0',
            'mom2': '0', 'proband2': 'mom2', 'dad2': '0', 'sib': 'mom2'})
        self.assertEqual(fathers, {'mom': '0', 'proband': 'dad', 'dad': '0',
            'mom2': '0', 'proband2': 'dad2', 'dad2': '0', 'sib': 'dad2'})
        self.assertEqual(children, {'proband': 'fam_ID',
            "proband2": "fam_ID2", "sib": "fam_ID2"})
        self.assertEqual(affected, {'mom': '1', 'proband': '2', 'dad': '1',
            'mom2': '1', 'proband2': '2', 'dad2': '1', 'sib': '2'})
        self.assertEqual(sex, {'mom': 'F', 'proband': 'F', 'dad': 'M',
            'mom2': 'F', 'proband2': 'F', 'dad2': 'M', 'sib': 'F'})
        self.assertEqual(vcfs, {'mom': '/path/to/mom_vcf.gz',
            'proband': '/path/to/proband_vcf.gz',
            'dad': '/path/to/dad_vcf.gz',
            'proband2': '/path/to/proband2_vcf.gz',
            'dad2': '/path/to/dad2_vcf.gz',
            'mom2': '/path/to/mom2_vcf.gz',
            'sib': '/path/to/sib_vcf.gz'})
    
    def test_load_families(self):
        """ check that load_families works correctly
        """
        
        # construct a temporary family that will have the same sample IDs etc
        # as for the one loaded from the ped file.
        family = Family("fam_ID")
        family.add_child("proband", "/path/to/proband_vcf.gz", "2", "F")
        family.add_mother("mom", "/path/to/mom_vcf.gz", "1", "F")
        family.add_father("dad", "/path/to/dad_vcf.gz", "1", "M")
        
        # load the ped file, and check that the load_families function returns
        # the expected Family object
        families = load_families(self.path)
        self.assertEqual(families, {"fam_ID": family})
        
        # add an extra family, with multiple sibs
        self.tempfile.write("fam_ID2  proband2 dad2  mom2  F  2  /path/to/proband2_vcf.gz\n")
        self.tempfile.write("fam_ID2  dad2     0     0     M  1  /path/to/dad2_vcf.gz\n")
        self.tempfile.write("fam_ID2  mom2     0     0     F  1  /path/to/mom2_vcf.gz\n")
        self.tempfile.write("fam_ID2  sib      dad2  mom2  F  2  /path/to/sib_vcf.gz\n")
        self.tempfile.flush()
        
        # construct a temporary family that will have the same sample IDs etc
        # as for the one loaded from the ped file.
        fam2 = Family("fam_ID2")
        fam2.add_child("proband2", "/path/to/proband2_vcf.gz", "2", "F")
        fam2.add_child("sib", "/path/to/sib_vcf.gz", "2", "F")
        fam2.add_mother("mom2", "/path/to/mom2_vcf.gz", "1", "F")
        fam2.add_father("dad2", "/path/to/dad2_vcf.gz", "1", "M")
        
        # load the ped file, and check that the load_families function returns
        # the expected Families objects
        families = load_families(self.path)
        self.assertEqual(set(families.values()), set([family, fam2]))
    
