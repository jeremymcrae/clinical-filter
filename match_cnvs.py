""" class to match CNVs based on their overlap
"""


class MatchCNVs(object):
    """ class to find if a parents CNV matches any of a childs CNVs
    """
    
    def __init__(self, variants):
        """ initiate the class with a dict of variants
        """
        
        self.variants = variants
        
        # figure out which of the childs variants are CNVs
        self.cnvs = []
        for key in self.variants:
            if len(key) == 3: # ignore SNVs, which are only (chrom, position)
                self.cnvs.append(key)
    
    def has_match(self, var):
        """ checks if a CNV has any matches amongst the childs CNVs
        
        Args:
            var: CNV object
        
        Returns:
            returns true if any of the childs CNVs are good matches
        """
        
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
    
    def get_overlap_key(self, var_key):
        """ returns the tuple for an overlapping CNV
        """
        
        if var_key in self.overlap:
            return self.overlap[var_key]
        else:
            return None



