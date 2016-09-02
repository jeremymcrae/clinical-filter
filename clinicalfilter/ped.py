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

import os
import sys

class Person(object):
    """creates an object for a person, with their ID, and VCF path
    """
    
    male_codes = set(["1", "m", "M", "male"])
    female_codes = set(["2", "f", "F", "female"])
    
    def __init__(self, family_id, person_id, dad_id, mom_id, sex, status, path):
        self.family_id = family_id
        self.person_id = person_id
        self.mom_id = mom_id
        self.dad_id = dad_id
        self.vcf_path = path
        self.sex = sex
        self.status = status
        
        # set a flag so we can check whether the child has been analysed
        self.analysed = False
    
    def __repr__(self):
        return 'Person(family_id="{}", person_id="{}", dad_id="{}", ' \
            'mom_id="{}", sex="{}", status="{}", path="{}")'.format(self.family_id,
            self.get_id(), self.dad_id, self.mom_id, self.get_gender(),
             self.get_affected_status(), self.get_path())
    
    def get_id(self):
        """returns the ID for a person.
        """
        return self.person_id
    
    def get_path(self):
        """returns the path to the VCF file for a person.
        """
        return self.vcf_path
    
    def get_affected_status(self):
        """returns the affected status for a person as a string
        """
        return self.status
    
    def is_affected(self):
        """returns true or false for affected, rather than the string value
        """
        # change how the affected status is encoded. Current DDD ped files
        # encode "1" for unaffected, and "2" for affected. Change this to
        # True/False values, and catch any unknown affected statuses.
        if self.status not in set(["1", "2"]):
            raise ValueError("unknown status: " + self.status + ", \
                should be 1: unaffected, 2: affected")
        
        return self.status == "2"
    
    def set_analysed(self):
        """ sets an individual as having been analysed
        """
        self.analysed = True
    
    def is_analysed(self):
        """ checks whether the individual has been analysed
        """
        return self.analysed
    
    def get_gender(self):
        """returns the gender for a person (1, M = male, 2, F = female).
        """
        return self.sex
    
    def is_male(self):
        """ returns True/False for whether the person is male
        """
        
        return self.get_gender() in self.male_codes
    
    def is_female(self):
        """ returns True/False for whether the person is male
        """
        
        return self.get_gender() in self.female_codes
    
    def check_gender(self, gender_code):
        """ makes sure that the parents match their expected gender.
        
        Rather than returning true/false I've chosen to raise an error, as this
        is a problem with the input data that needs to be fixed.
        
        Args:
            gender_code: mothers, are "2", while fathers, are "1"
        """
        
        if self.is_male():
            current_gender_codes = self.male_codes
        elif self.is_female():
            current_gender_codes = self.female_codes
        else:
            raise ValueError("unknown gender code: " + self.get_gender())
        
        if gender_code not in current_gender_codes:
            raise ValueError(self.person_id + " is listed as gender " + \
                self.get_gender() + ", which differs from the sex expected " + \
                "as a parent)")
    
    def __gt__(self, other):
        """ implement greater than check, for sorting children in a Family
        """
        
        return self.get_id() > other.get_id()
    
    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
    
    def __hash__(self):
        """ get a unique hash for the object from the sample strings
        """
        
        return hash(tuple([self.get_id(), self.get_path(),
            self.get_affected_status(), self.get_gender()]))

class Family(object):
    """creates a family, with VCF paths, IDs, and affected statuses
    """
    
    def __init__(self, family_id, children=None, mother=None, father=None):
        """ initiates the class with the ID for the family
        """
        self.family_id = family_id
        self.children = children
        if self.children is None:
            self.children = []
        
        self.father = father
        self.mother = mother
        self.set_child()
    
    def __repr__(self):
        return 'Family(family_id="{}", children={}, mother={}, ' \
            'father={})'.format(self.family_id, self.children, self.mother,
            self.father)
    
    def __iter__(self):
        for member in self.children + [self.mother, self.father]:
            yield member
    
    def __gt__(self, other):
        return self.family_id > other.family_id
    
    def has_parents(self):
        """ returns True/False for whether the family includes parental info
        
        Currently requires both parents to be included
        """
        
        return self.father is not None and self.mother is not None
    
    def add_child(self, sample_id, dad_id, mom_id, sex, status, path):
        """ adds a child
        
        Args:
            sample_id: individual ID string
            path: path to childs VCF file
            affected_status: affected status string for child
            gender: gender string for child
        """
        child = Person(self.family_id, sample_id, dad_id, mom_id, sex, status, path)
        self.children.append(child)
    
    def add_mother(self, sample_id, mom_id, dad_id, sex, status, path):
        # raise an error if we try to add a different mother to the family
        if self.mother is not None:
            if sample_id != self.mother.get_id():
                raise ValueError(self.family_id, "already has a mother")
        
        self.mother = Person(self.family_id, sample_id, dad_id, mom_id, sex, status, path)
        self.mother.check_gender("2")
    
    def add_father(self, sample_id, mom_id, dad_id, sex, status, path):
        # raise an error if we try to add a different father to the family
        if self.father is not None:
            if sample_id != self.father.get_id():
                raise ValueError(self.family_id, "already has a father")
        
        self.father = Person(self.family_id, sample_id, dad_id, mom_id, sex, status, path)
        self.father.check_gender("1")
    
    def set_child(self):
        """ define the child to be examined
        """
        for child in self.children:
            if not child.is_analysed():
                self.child = child
                return
        
        # if we have run through all the children, set the child to None
        self.child = None
    
    def set_child_examined(self):
        """ once a child has been examined, mark it as such in the children list
        """
        for child_position in range(len(self.children)):
            child = self.children[child_position]
            if child.get_id() == self.child.get_id():
                self.children[child_position].set_analysed()
        
        self.set_child()
    
    def __eq__(self, other):
        """ check for equality between two Family objects
        """
    
        return hash(self) == hash(other)
    
    def __hash__(self):
        """ construct a unique hash, based on the hashes for the family members
        """
        
        parts = tuple([self.family_id, hash(self.mother), hash(self.father)] + \
            [hash(x) for x in sorted(self.children)])
        
        return hash(parts)

def open_ped(path):
    """ opens a ped file, and groups individuals into families
    
    The PED file is in LINKAGE PED format, with the first six columns
    specfifying the individual and how they are related to other individuals. In
    contrast to other PED files, the genotypes are specified as a path to a VCF
    file for the individual.
    
    Args:
        path: path to the ped file
    
    Returns:
        dictionary of lists of Person objects for individuals per family,
        indexed by family ID.
    """
    
    if not os.path.exists(path):
        sys.exit("Path to ped file does not exist: " + path)
    
    # group the lines in the family relationships file by family
    families = {}
    with open(path) as handle:
        for line in handle:
            # parse the line as a Person object, to assist downstream organising
            line = Person(*line.strip().split())
            fam_id = line.family_id
            
            # add the Person to a list of lines for the family
            if fam_id not in families:
                families[fam_id] = []
            
            families[fam_id].append(line)
    
    return families

def load_families(path):
    """ Creates a list of family data from a PED file.
    
    Args:
        path: path to the ped file
        
    Returns:
        list of Family objects
    """
    
    family_lines = open_ped(path)
    
    families = []
    for fam_id, lines in family_lines.items():
        children = lines
        parents = []
        if len(lines) > 1:
            children = [ x for x in lines if x.dad_id != '0' or x.mom_id != '0' ]
            parents = [ x for x in lines if x.dad_id == '0' or x.mom_id == '0' ]
        
        moms = [ x for x in parents if x.get_id() == children[0].mom_id ]
        dads = [ x for x in parents if x.get_id() == children[0].dad_id ]
        
        mom = None
        if len(moms) == 1:
            mom = moms[0]
        
        dad = None
        if len(dads) == 1:
            dad = dads[0]
        
        family = Family(fam_id, children=children, mother=mom, father=dad)
        families.append(family)
    
    return families
