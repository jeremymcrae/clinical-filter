""" class to match CNVs based on their overlap
"""

import math

class MatchCNVs(object):
    """ class to find if a CNV matches any of another individuals CNVs
    """
    
    def __init__(self, variants):
        """ initiate the class with a dict of variants
        """
        
        # figure out which of the variants are CNVs
        self.cnvs = []
        for key in variants:
            if len(key) == 3: # ignore SNVs, which are only (chrom, position)
                self.cnvs.append(key)
        
    def has_match(self, var):
        """ checks if any of the individuals CNVs overlap the current CNV
        
        Args:
            var: CNV object
        
        Returns:
            returns true if any of the individuals CNVs overlap
        """
        
        var_key = var.get_key() # CNVs are (chrom, start, end)
        var_chrom = var_key[0]
        var_start = int(var_key[1])
        var_end = int(var_key[2])
        
        self.overlap = {}
        for key in self.cnvs:
            
            chrom = key[0]
            start = key[1]
            end = key[2]
            
            if var_chrom != chrom:
                continue
            
            # check if the individuals CNVs end points lie within the CNV that
            # we are examining, or if the individuals CNV surrounds the examined
            # CNV
            if var_start <= int(end) and var_end >= int(start) and self.similar_size(var, start, end):
                self.overlap[var_key] = key
        
        return len(self.overlap) > 0
    
    def calculate_cnv_size_tolerance(self, var):
        """ calculates the size range of CNVs that might match a given CNV size.
        
        Args:
            var: CNV object
        
        Returns:
            tuple of minimum and maximum sizes in base pairs
        """
        
        var_key = var.get_key()
        var_start = int(var_key[1])
        var_end = int(var_key[2])
        size = var_end - var_start
        
        min_size = size - abs(100 * math.sqrt(size + 2500)) + 5000
        max_size = size + 100 * math.sqrt(size)
        
        return (min_size, max_size)
    
    def similar_size(self, var, start, end):
        """ checks if the current CNV matches the overlapping child CNVs size
        
        Args:
            var: CNV object
        
        Returns:
            returns true if any of the overlapping CNVs have similar sizes
        """
        
        (min_size, max_size) = self.calculate_cnv_size_tolerance(var)
        
        size = int(end) - int(start)
        
        if max_size > size > min_size:
            return True
        
        return False
    
    def get_overlap_key(self, var_key):
        """ returns the tuple for an overlapping CNV
        """
        
        if var_key in self.overlap:
            return self.overlap[var_key]
        else:
            return None



