""" class for filtering Exome based CNV calls
"""

class BreakdancerCNV(object):
    """ class to filter breakdancer genome CNV calls
    """
    
    def __init__(self, cnv):
        """ initialise the class with a CNV
        """
        self.cnv = cnv
    
    def filter_cnv(self, track_variant):
        """ filters the CNV
        """
        
        return True
        
