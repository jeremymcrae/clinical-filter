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
from clinicalfilter.ped import Person


class TestPerson(unittest.TestCase):
    """
    """
    
    def setUp(self):
        """ define a default Person object
        """
        
        family_id = "fam_id"
        person_id = "name"
        path = "/home/filename.vcf"
        status = "2"
        sex = "1"
        mom_id = "mom_id"
        dad_id = "dad_id"
        
        self.person = Person(family_id, person_id, dad_id, mom_id, sex, status, path)
    
    def test_get_id(self):
        """ test that get_id() works correctly
        """
        
        self.person.person_id = "test_id"
        self.assertEqual(self.person.get_id(), "test_id")
    
    def test_get_path(self):
        """ test that get_path() works correctly
        """
        
        self.person.vcf_path = "test_path"
        self.assertEqual(self.person.get_path(), "test_path")
    
    def test_get_affected_status(self):
        """ test that get_affected_status works correctly
        """
        
        self.person.status = "test_status"
        self.assertEqual(self.person.get_affected_status(), "test_status")
    
    def test_is_affected(self):
        """ test that is_affected() works correctly
        """
        
        # check that a status of "1" means unaffected
        self.person.status = "1"
        self.assertFalse(self.person.is_affected())
        
        # check that a status of "2" means affected
        self.person.status = "2"
        self.assertTrue(self.person.is_affected())
        
        # check that statuses other than "1" or "2" raise an error
        self.person.status = "3"
        with self.assertRaises(ValueError):
            self.assertFalse(self.person.is_affected())
    
    def test_set_analysed(self):
        """ test that set_analysed() and is_analysed() work correctly
        """
        
        # check that by default a person is not analysed
        self.assertFalse(self.person.is_analysed())
        
        # check that set_analysed() changes the flag correctly
        self.person.set_analysed()
        self.assertTrue(self.person.is_analysed())
        
        # check that repeating set_analysed() does not affected anything
        self.person.set_analysed()
        self.assertTrue(self.person.is_analysed())
    
    def test_get_gender(self):
        """ test that get_gender() works correctly
        """
        
        self.person.sex = "M"
        self.assertEqual(self.person.get_gender(), "M")
    
    def test_is_male(self):
        """ test that is_male() works correctly
        """
        
        male_codes = ["1", "M", "m", "male"]
        for code in male_codes:
            self.person.sex = code
            self.assertTrue(self.person.is_male())
        
        female_codes = ["2", "f", "F", "female"]
        for code in female_codes:
            self.person.sex = code
            self.assertFalse(self.person.is_male())
    
    def test_is_female(self):
        """ test that is_female() works correctly
        """
        
        female_codes = ["2", "f", "F", "female"]
        for code in female_codes:
            self.person.sex = code
            self.assertTrue(self.person.is_female())
        
        male_codes = ["1", "M", "m", "male"]
        for code in male_codes:
            self.person.sex = code
            self.assertFalse(self.person.is_female())
    
    def test_check_gender(self):
        """ test that check_gender() works correctly
        """
        
        # check that a male gender doesn't raise an error if we check that it
        # is male
        self.person.sex = "M"
        self.person.check_gender("1")
        
        # check that we raise an error if the gender doesn't match the expected
        # (useful in catching errors with the parents)
        with self.assertRaises(ValueError):
            self.person.check_gender("2")
        
        # check that a female gender doesn't raise an error if we check that it
        # is female
        self.person.sex = "F"
        self.person.check_gender("2")
        
        # check that we raise an error if the gender doesn't match the expected
        # (useful in catching errors with the parents)
        with self.assertRaises(ValueError):
            self.person.check_gender("1")
        
        # check that we raise an error for nonstandard gender codes
        self.person.sex = "NA"
        with self.assertRaises(ValueError):
            self.person.check_gender("2")


if __name__ == '__main__':
    unittest.main()
