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

from clinicalfilter.ped import open_ped, load_families, Person, Family

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
    
    def test_open_ped_single_family(self):
        """ check that we correctly parse a ped file with a single trio
        """
        
        # load all the components from the file
        self.assertEqual(open_ped(self.path), {'fam_ID': [
            Person(family_id="fam_ID", person_id="proband", dad_id="dad",
                mom_id="mom", sex="F", status="2", path="/path/to/proband_vcf.gz"),
            Person(family_id="fam_ID", person_id="dad", dad_id="0", mom_id="0",
                sex="M", status="1", path="/path/to/dad_vcf.gz"),
            Person(family_id="fam_ID", person_id="mom", dad_id="0", mom_id="0",
                sex="F", status="1", path="/path/to/mom_vcf.gz")]})
    
    def test_open_ped_multiple_sibs(self):
        """ check that we correctly parse a ped file with multiple siblings
        """
        
        # add an extra sibling
        self.tempfile.write("fam_ID  sib   dad  mom  F  2  /path/to/sib_vcf.gz\n")
        self.tempfile.flush()
        
        # load all the components from the file
        self.assertEqual(open_ped(self.path), {'fam_ID': [
            Person(family_id="fam_ID", person_id="proband", dad_id="dad",
                mom_id="mom", sex="F", status="2", path="/path/to/proband_vcf.gz"),
            Person(family_id="fam_ID", person_id="dad", dad_id="0", mom_id="0",
                sex="M", status="1", path="/path/to/dad_vcf.gz"),
            Person(family_id="fam_ID", person_id="mom", dad_id="0", mom_id="0",
                sex="F", status="1", path="/path/to/mom_vcf.gz"),
            Person(family_id="fam_ID", person_id="sib", dad_id="mom", mom_id="dad",
                sex="F", status="2", path="/path/to/sib_vcf.gz")]})
    
    def test_open_ped_multiple_families(self):
        """ check that we correctly parse a ped file with multiple families
        """
        
        # add an extra family, with multiple sibs
        self.tempfile.write("fam_ID2  proband2  dad2  mom2  F  2  /path/to/proband2_vcf.gz\n")
        self.tempfile.write("fam_ID2  dad2      0     0     M  1  /path/to/dad2_vcf.gz\n")
        self.tempfile.write("fam_ID2  mom2      0     0     F  1  /path/to/mom2_vcf.gz\n")
        self.tempfile.write("fam_ID2  sib       dad2  mom2  F  2  /path/to/sib_vcf.gz\n")
        self.tempfile.flush()
        
        # load all the components from the file
        self.assertEqual(open_ped(self.path), {
            'fam_ID': [
                Person(family_id="fam_ID", person_id="proband", dad_id="dad",
                    mom_id="mom", sex="F", status="2", path="/path/to/proband_vcf.gz"),
                Person(family_id="fam_ID", person_id="dad", dad_id="0", mom_id="0",
                    sex="M", status="1", path="/path/to/dad_vcf.gz"),
                Person(family_id="fam_ID", person_id="mom", dad_id="0", mom_id="0",
                    sex="F", status="1", path="/path/to/mom_vcf.gz")],
            'fam_ID2': [
                Person(family_id="fam_ID2", person_id="proband2", dad_id="dad2",
                    mom_id="mom2", sex="F", status="2", path="/path/to/proband2_vcf.gz"),
                Person(family_id="fam_ID2", person_id="dad2", dad_id="0", mom_id="0",
                    sex="M", status="1", path="/path/to/dad2_vcf.gz"),
                Person(family_id="fam_ID2", person_id="mom2", dad_id="0", mom_id="0",
                    sex="F", status="1", path="/path/to/mom2_vcf.gz"),
                Person(family_id="fam_ID2", person_id="sib", dad_id="dad2",
                    mom_id="mom2", sex="F", status="2", path="/path/to/sib_vcf.gz"),
            ]})
    
    def test_load_families(self):
        """ check that load_families works correctly
        """
        
        # construct a temporary family that will have the same sample IDs etc
        # as for the one loaded from the ped file.
        family = Family("fam_ID")
        family.add_child("proband", 'dad', 'mom', 'F', '2', "/path/to/proband_vcf.gz")
        family.add_mother("mom", '0', '0', 'F', '1', "/path/to/mom_vcf.gz")
        family.add_father("dad", '0', '0', 'M', '1',  "/path/to/dad_vcf.gz")
        
        # load the ped file, and check that the load_families function returns
        # the expected Family object
        self.assertEqual(load_families(self.path), [family])
        
        # add an extra family, with multiple sibs
        self.tempfile.write("fam_ID2  proband2 dad2  mom2  F  2  /path/to/proband2_vcf.gz\n")
        self.tempfile.write("fam_ID2  dad2     0     0     M  1  /path/to/dad2_vcf.gz\n")
        self.tempfile.write("fam_ID2  mom2     0     0     F  1  /path/to/mom2_vcf.gz\n")
        self.tempfile.write("fam_ID2  sib      dad2  mom2  F  2  /path/to/sib_vcf.gz\n")
        self.tempfile.flush()
        
        # construct a temporary family that will have the same sample IDs etc
        # as for the one loaded from the ped file.
        fam2 = Family("fam_ID2")
        fam2.add_child("proband2", 'dad2', 'mom2', 'F', '2', "/path/to/proband2_vcf.gz")
        fam2.add_child("sib", 'dad2', 'mom2', 'F', '2', "/path/to/sib_vcf.gz")
        fam2.add_mother("mom2", '0', '0', 'F', '1', "/path/to/mom2_vcf.gz")
        fam2.add_father("dad2", '0', '0', 'M', '1', "/path/to/dad2_vcf.gz")
        
        # load the ped file, and check that the load_families function returns
        # the expected Families objects
        self.assertEqual(sorted(load_families(self.path)), sorted([family, fam2]))
    
