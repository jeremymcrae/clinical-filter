""" A class for loading ped files, and sorting the lines into families
"""

import os
import sys

class Person(object):
    """creates an object for a person, with their ID, and VCF path
    """
    
    male_codes = set(["1", "m", "M", "male"])
    female_codes = set(["2", "f", "F", "female"])
    
    def __init__(self, person_id, vcf_path, affected_status, gender):
        self.person_id = person_id
        self.vcf_path = vcf_path
        self.gender = gender
        self.affected_status = affected_status
        
        # set a flag so we can check whether the child has been analysed
        self.analysed = False
    
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
        return self.affected_status
    
    def is_affected(self):
        """returns true or false for affected, rather than the string value
        """
        # change how the affected status is encoded. Current DDD ped files
        # encode "1" for unaffected, and "2" for affected. Change this to
        # True/False values, and catch any unknown affected statuses.
        if self.affected_status not in set(["1", "2"]):
            raise ValueError("unknown status: " + self.affected_status + ", \
                should be 1: unaffected, 2: affected")
        
        return self.affected_status == "2"
    
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
        return self.gender
    
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
    
    def __hash__(self):
        """ get a unique hash for the object from the sample strings
        """
        
        return hash(tuple([self.get_id(), self.get_path(),
            self.get_affected_status(), self.get_gender()]))

class Family(object):
    """creates a family, with VCF paths, IDs, and affected statuses
    """
    
    def __init__(self, family_id):
        """ initiates the class with the ID for the family
        """
        self.family_id = family_id
        self.children = []
        self.father = None
        self.mother = None
        self.child = None
    
    def has_parents(self):
        """ returns True/False for whether the family includes parental info
        
        Currently requires both parents to be included
        """
        
        return self.father is not None and self.mother is not None
    
    def add_child(self, sample_id, path, affected_status, gender):
        """ adds a child
        
        Args:
            sample_id: individual ID string
            path: path to childs VCF file
            affected_status: affected status string for child
            gender: gender string for child
        """
        child = Person(sample_id, path, affected_status, gender)
        self.children.append(child)
    
    def add_mother(self, sample_id, path, affected_status, gender):
        # raise an error if we try to add a different mother to the family
        if self.mother is not None:
            if sample_id != self.mother.get_id():
                raise ValueError(self.family_id, "already has a mother")
        
        self.mother = Person(sample_id, path, affected_status, gender)
        self.mother.check_gender("2")
    
    def add_father(self, sample_id, path, affected_status, gender):
        # raise an error if we try to add a different father to the family
        if self.father is not None:
            if sample_id != self.father.get_id():
                raise ValueError(self.family_id, "already has a father")
        
        self.father = Person(sample_id, path, affected_status, gender)
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

def load_ped_file(path):
    """Loads a PED file containing details for multiple trios.
    
    The PED file is in LINKAGE PED format, with the first six columns
    specfifying the indivudual and how they are related to toher individuals. In
    contrast to other PED files, the genotypes are specified as a path to a VCF
    file for the individual.
    
    Args:
        path: path to the ped file
    
    Returns:
        mothers: dictionary of maternal IDs, indexed by the childs ID
        fathers: dictionary of paternal IDs, indexed by the childs ID
        children: dictionary of family IDs, indexed by the childs ID
        affected: dictionary of affected statuses, indexed by individual ID
        sex: dictionary of genders, indexed by individual ID
        vcfs: dictionary of VCF paths, indexed by individual ID
    """
    
    if not os.path.exists(path):
        sys.exit("Path to ped file does not exist: " + path)
    
    mothers = {}
    fathers = {}
    children = {}
    affected = {}
    sex = {}
    vcfs = {}
    
    ped = open(path, "r")
    for line in ped:
        line = line.strip().split()
        
        family_id = line[0]
        individual_id = line[1]
        paternal_id = line[2]
        maternal_id = line[3]
        gender = line[4]
        affected_status = line[5]
        path = line[6]
        
        # make sure we can match individuals to their paths, and affected status
        vcfs[individual_id] = path
        affected[individual_id] = affected_status
        sex[individual_id] = gender
        
        # track the child, maternal and paternal IDs
        if paternal_id != "0" or maternal_id != "0":
            children[individual_id] = family_id
        
        fathers[individual_id] = paternal_id
        mothers[individual_id] = maternal_id
    
    ped.close()
    
    return mothers, fathers, children, affected, sex, vcfs

def load_families(path):
    """Creates a dictionary of family data from a PED file.
    
    Args:
        path: path to the ped file
        
    Returns:
        families: a dictionary of Family objects, indexed by family IDs
    """
    
    # do an initial parse of the ped file
    mothers, fathers, children, affected, sex, vcfs = load_ped_file(path)
    
    families = {}
    
    # put all the family info into a trio class
    for child_id, family_id in children.items():
        # if the family hasn't been included already, generate a new trio
        if family_id not in families:
            families[family_id] = Family(family_id)
        
        trio = families[family_id]
        father = fathers[child_id]
        mother = mothers[child_id]
        
        # add the child to the family, and set the child to be examined
        trio.add_child(child_id, vcfs[child_id], affected[child_id], sex[child_id])
        trio.set_child()
        
        # add parents, but allow for children without parents listed in the ped
        # file
        if father in vcfs and father in affected:
            trio.add_father(father, vcfs[father], affected[father], sex[father])
        if mother in vcfs and mother in affected:
            trio.add_mother(mother, vcfs[mother], affected[mother], sex[mother])
        
        families[family_id] = trio
    
    return families
