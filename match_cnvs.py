""" class to match CNVs based on their overlap and relative sizes and positions
"""


class MatchCNVs(object):
    """ class to find if a parents CNV matches any of a childs CNVs
    """
    
    def __init__(self, child_variants):
        """ initiate the class with the childs variants
        """
        
        self.child_variants = child_variants
        
        # figure out which of the childs variants are CNVs
        self.cnvs = []
        for key in self.child_variants:
            if len(key) == 3: # ignore SNVs, which are only (chrom, position)
                self.cnvs.append(key)
    
    def has_match(self, var):
        """ checks if a CNV has any matches amongst the childs CNVs
        
        Args:
            var: CNV object
        
        Returns:
            returns true if any of the childs CNVs are good matches
        """
        
        # TODO: swap to code that compares CNV sizes in line below
        # if self.any_overlap(var) and self.similar_size(var):
        if self.any_overlap(var):
            return True
        else:
            return False
        
    def any_overlap(self, var):
        """ checks if any of the childs CNVs overlap the current CNV
        
        Args:
            var: CNV object
        
        Returns:
            returns true if any of the childs CNVs overlap
        """
        
        var_key = var.get_key() # CNVs are (chrom, start, end)
        var_chrom = var_key[0]
        var_start = int(var_key[1])
        var_end = int(var_key[2])
        
        self.overlap = {}
        for key in self.cnvs:
            
            child_chrom = key[0]
            child_start = int(key[1])
            child_end = int(key[2])
            
            if var_chrom != child_chrom:
                continue
            
            # check if the childs CNVs end points lie within the CNV that we
            # are examining, or if the childs CNV surrounds the examined CNV
            if var_end >= child_start >= var_start or \
               var_end >= child_end >= var_start or \
               child_end >= var_start >= child_start:
                self.overlap[var_key] = key
        
        if len(self.overlap) > 0:
            return True
        else:
            return False
        
    def similar_size(self, var):
        """ checks if the current CNV matches the overlapping child CNVs size
        
        Args:
            var: CNV object
        
        Returns:
            returns true if any of the overlapping CNVs have similar sizes
        """
        (min_size, max_size) = var.calculate_cnv_size_tolerance()
        
        for overlap in self.overlap:
            var_chrom = overlap[0]
            var_start = int(overlap[1])
            var_end = int(overlap[2])
            var_size = var_end - var_start
            
            if max_size > var_size > min_size:
                self.overlap = overlap
                return True
        
        return False
    
    def get_overlap_key(self, var_key):
        """ returns the tuple for an overlapping CNV
        """
        
        if var_key in self.overlap:
            return self.overlap[var_key]
        else:
            return None



