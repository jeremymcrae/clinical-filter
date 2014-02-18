""" class for filtering Exome based CNV calls
"""

class ExomeCNV(object):
    """ class to filter exome CNV calls
    """
    
    def __init__(self, cnv):
        """ initialise the class with a CNV
        """
        self.cnv = cnv
    
    def filter_cnv(self, track_variant):
        """ filters the CNV
        """
        
        passes = True
        if self.fails_convex_score():
            passes = False
            if track_variant:
                print("failed CONVEX score", self.cnv.info["CONVEX"])
        elif self.fails_population_frequency():
            passes = False
            if track_variant:
                print("failed pop freq", self.cnv.info["RC50INTERNALFREQ"])
        elif self.fails_mad_ratio():
            passes = False
            if track_variant:
                print("failed mad ratio", self.cnv.info["MEANLR2"], self.cnv.info["MADL2R"])
        elif self.fails_commmon_forwards():
            passes = False
            if track_variant:
                print("failed commonforwards", self.cnv.info["COMMONFORWARDS"])
        # elif self.fails_no_exons():
        #     passes = False
        #     if track_variant:
        #         print("failed no exons", self.info["NUMBEREXONS"])
        
        return passes
    
    def fails_convex_score(self):
        """ checks if the convex score is out of bounds
        """
        return float(self.cnv.info["CONVEXSCORE"]) <= 7
    
    def fails_population_frequency(self):
        """ checks if the population frequency for the CNV is too high
        """
        
        return float(self.cnv.info["RC50INTERNALFREQ"]) > 0.01
    
    def fails_mad_ratio(self):
        """ checks if the MAD ratio is too low
        """
        
        try:
            return abs(float(self.cnv.info["MEANLR2"])/float(self.cnv.info["MADL2R"])) < 5
        except ZeroDivisionError:
            return True
    
    def fails_commmon_forwards(self):
        """ checks if the COMMONFORWARDS value is too high
        """
        
        return float(self.cnv.info["COMMONFORWARDS"]) > 0.8
    
    def fails_no_exons(self):
        """ checks that the CNV overlaps at least one exon
        """
        
        return float(self.cnv.info["NUMBEREXONS"]) < 1
