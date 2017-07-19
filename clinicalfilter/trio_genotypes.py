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

class TrioGenotypes(object):
    """ loads variant data from individuals into a single object, so we can easily
    get genotype data for a single variant from all the family members.
    """
    
    def __init__(self, chrom=None, pos=None, child=None, mother=None,
            father=None, debug_chrom=None, debug_pos=None):
        """ initiate the class with the childs variant
        
        Args:
            child_variant: Variant object
        """
        
        self.chrom = chrom
        self.pos = pos
        
        self.child = child
        self.mother = mother
        self.father = father
        
        self.debug_chrom = debug_chrom
        self.debug_pos = debug_pos
    
    def get_chrom(self):
        if self.child is not None:
            return self.child.get_chrom()
        else:
            return self.chrom
    
    def get_position(self):
        if self.child is not None:
            return self.child.get_position()
        else:
            return self.pos
    
    def get_genes(self):
        if self.child is not None:
            return self.child.info.get_genes()
        else:
            return None
    
    def get_range(self):
        if self.child is not None:
            return self.child.info.get_range()
        else:
            return None
    
    def is_cnv(self):
        if self.child is not None:
            return self.child.is_cnv()
        else:
            return None
    
    def chrom_to_int(self, chrom):
        """ converts a chromosome string to an int (if possible) for sorting.
        
        Args:
            chrom: string (eg "1", "2", "3" ... "22", "X", "Y")
            
        Returns:
            int value of chrom
        """
        
        # set the integer values of the sex chromosomes
        chrom_dict = {"X": 23, "CHRX": 23, "Y": 24, "CHRY": 24, "MT": 25, \
            "CHRMT": 25}
            
        try:
            chrom = int(chrom)
        except ValueError:
            chrom = chrom_dict[chrom.upper()]
        
        return chrom
    
    def __eq__(self, other):
        return self.__hash__() == other.__hash__()
    
    def __lt__(self, other):
        return (self.chrom_to_int(self.get_chrom()), int(self.get_position())) < \
              (self.chrom_to_int(other.get_chrom()), int(other.get_position()))
    
    def __repr__(self):
        return 'TrioGenotypes(chrom="{}", pos={}, child={}, mother={},' \
            'father={})'.format(self.get_chrom(), self.get_position(),
            self.child, self.mother, self.father)
    
    def __str__(self):
        genotype = '{}/{}/{}'.format(*self.get_trio_genotype())
        genotype = genotype.replace('None', 'NA')
        
        return '{}:{} - {}'.format(self.get_chrom(), self.get_position(), genotype)
    
    def __hash__(self):
        return hash(self.__repr__())
    
    def get_inheritance_type(self):
        if self.child is not None:
            return self.child.get_inheritance_type()
        else:
            return None
    
    def get_trio_genotype(self):
        child = self.child.get_genotype()
        
        mother, father = None, None
        if self.mother is not None:
            mother = self.mother.get_genotype()
        
        if self.father is not None:
            father = self.father.get_genotype()
        
        return (child, mother, father)
    
    def get_de_novo_genotype(self):
        """ get the de novo genotype combination for the chromosome/sex
        """
        
        # set the standard de novo genotype combination
        de_novo_genotype = (1, 0, 0)
        
        # account for X chrom de novos in males
        if self.get_inheritance_type() == "XChrMale":
            de_novo_genotype = (2, 0, 0)
        
        return de_novo_genotype
    
    def passes_de_novo_checks(self, pp_filter):
        """ checks if the child's de novo variants passes filters
        
        Some variants are de novo in the child, and if they are, then we should
        subject them to additional filtering to see if they have been passed by
        denovogear. Note that both parents have genotypes specified, as we
        have checked at an earlier stage that the parents are present.
        
        Args:
            pp_filter: float between 0 and 1, being the threshold for the PP_DNM filter
        
        Returns:
            boolean value for whether the variant should be included
        """
        
        # currently hard code the filtering fields. The de novo field indicates
        # whether the variant is de novo, I don't know how this is assigned. The
        # project filter field indicates an internal filter, curently whether
        # the variant passed MAF, alternate frequency, and segmental duplication
        # criteria.
        de_novo_field = set(["DENOVO-SNP", "DENOVO-INDEL"])
        
        # if the variant is not de novo, don't worry about de novo filtering
        if self.get_trio_genotype() != self.get_de_novo_genotype():
            return True
        
        # check the VCF record to see whether the variant has been screened out.
        # Either DENOVO-SNP or DENOVO-INDEL should be in the info.
        if len(set(self.child.info) & de_novo_field) < 1:
            if self.get_chrom() == self.debug_chrom and self.get_position() == self.debug_pos:
                print(self, "failed DENOVO-SNP/INDEL check")
            return False
        
        if "PP_DNM" in self.child.format and \
                float(self.child.format["PP_DNM"]) < pp_filter:
            if self.get_chrom() == self.debug_chrom and self.get_position() == self.debug_pos:
                print(self, "failed PP_DNM threshold")
            return False
        
        if "TEAM29_FILTER" in self.child.format:
            if self.child.format["TEAM29_FILTER"] != "PASS":
                if self.get_chrom() == self.debug_chrom and self.get_position() == self.debug_pos:
                    print(self, "failed TEAM29_FILTER")
                return False
        
        return True
