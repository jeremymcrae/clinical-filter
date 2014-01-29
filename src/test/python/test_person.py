""" unit testing of the Person class
"""

import unittest
from clinicalfilter.ped import Person


class TestPerson(unittest.TestCase):
    """
    """
    
    def setUp(self):
        """ define a default Person object
        """
        
        ID = "name"
        path = "/home/filename.vcf"
        status = "2"
        gender = "1"
        
        self.person = Person(ID, path, status, gender)
    
    def test_get_ID(self):
        """ test that get_ID() works correctly
        """
        
        self.person.person_ID = "test_id"
        self.assertEqual(self.person.get_ID(), "test_id")
    
    def test_get_path(self):
        """ test that get_path() works correctly
        """
        
        self.person.VCF_path = "test_path"
        self.assertEqual(self.person.get_path(), "test_path")
    
    def test_get_affected_status(self):
        """ test that get_affected_status works correctly
        """
        
        self.person.affected_status = "test_status"
        self.assertEqual(self.person.get_affected_status(), "test_status")
    
    def test_is_affected(self):
        """ test that is_affected() works correctly
        """
        
        # check that a status of "1" means unaffected
        self.person.affected_status = "1"
        self.assertFalse(self.person.is_affected())
        
        # check that a status of "2" means affected
        self.person.affected_status = "2"
        self.assertTrue(self.person.is_affected())
        
        # check that statuses other than "1" or "2" raise an error
        self.person.affected_status = "3"
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
        
        self.person.gender = "M"
        self.assertEqual(self.person.get_gender(), "M")
    
    def test_is_male(self):
        """ test that is_male() works correctly
        """
        
        male_codes = ["1", "M", "m", "male"]
        for code in male_codes:
            self.person.gender = code
            self.assertTrue(self.person.is_male())
        
        female_codes = ["2", "f", "F", "female"]
        for code in female_codes:
            self.person.gender = code
            self.assertFalse(self.person.is_male())
    
    def test_is_female(self):
        """ test that is_female() works correctly
        """
        
        female_codes = ["2", "f", "F", "female"]
        for code in female_codes:
            self.person.gender = code
            self.assertTrue(self.person.is_female())
        
        male_codes = ["1", "M", "m", "male"]
        for code in male_codes:
            self.person.gender = code
            self.assertFalse(self.person.is_female())
    
    def test_check_gender(self):
        """ test that check_gender() works correctly
        """
        
        # check that a male gender doesn't raise an error if we check that it 
        # is male
        self.person.gender = "M"
        self.person.check_gender("1")
        
        # check that we raise an error if the gender doesn't match the expected
        # (useful in catching errors with the parents)
        with self.assertRaises(ValueError):
            self.person.check_gender("2")
        
        # check that a female gender doesn't raise an error if we check that it 
        # is female
        self.person.gender = "F"
        self.person.check_gender("2")
        
        # check that we raise an error if the gender doesn't match the expected
        # (useful in catching errors with the parents)
        with self.assertRaises(ValueError):
            self.person.check_gender("1")
        
        # check that we raise an error for nonstandard gender codes
        self.person.gender = "NA"
        with self.assertRaises(ValueError):
            self.person.check_gender("2")
        
    
        
        
    
    




unittest.main()

