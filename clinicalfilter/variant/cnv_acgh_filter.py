""" class for filtering aCGH CNVs
"""

class ACGH_CNV(object):
    """ class for filtering array CGH CNV calls
    """
    
    def __init__(self, cnv):
        """ initialise the class with a CNV
        """
        
        self.cnv = cnv
    
    def filter_cnv(self, track_variant):
        """ filters the CNV
        """
        
        passes = True
        if self.fails_mad_ratio():
            passes = False
            if track_variant:
                print("failed mad ratio", self.cnv.info["MEANLR2"], self.cnv.info["MADL2R"])
        elif self.fails_wscore():
            passes = False
            if track_variant:
                print("failed wscore", self.cnv.info["WSCORE"])
        elif self.fails_callp():
            passes = False
            if track_variant:
                print("failed callp", self.cnv.info["CALLP"])
        elif self.fails_commmon_forwards():
            passes = False
            if track_variant:
                print("failed commonforwards", self.cnv.info["COMMONFORWARDS"])
        elif self.fails_meanlr2():
            passes = False
            if track_variant:
                print("failed meanlr2", self.cnv.info["MEANLR2"])
        elif self.fails_no_exons():
            passes = False
            if track_variant:
                print("failed no exons", self.cnv.info["NUMBEREXONS"])
        elif self.fails_frequency():
            passes = False
            if track_variant:
                print("failed frequency", self.cnv.info["ACGH_RC_FREQ50"])
        elif self.fails_cifer_inh():
            passes = False
            if track_variant:
                print("failed CIFER inheritance", self.cnv.format["CIFER_INHERITANCE"])
        
        return passes
    
    def fails_mad_ratio(self):
        """ checks if the MAD ratio is too low.
        
        Note that this filter has been lowered, so that it will no longer
        exclude any variants. This function could probably be removed.
        """
        
        try:
            return abs(float(self.cnv.info["MEANLR2"])/float(self.cnv.info["MADL2R"])) < 0
        except ZeroDivisionError:
            return True
        
    def fails_wscore(self):
        """ checks if the WSCORE value is too low
        """
        
        return float(self.cnv.info["WSCORE"]) < 0.4
    
    def fails_callp(self):
        """ checks if the CALLP value is too high
        """
        
        return float(self.cnv.info["CALLP"]) > 0.01
    
    def fails_commmon_forwards(self):
        """ checks if the COMMONFORWARDS value is too high
        """
        
        return float(self.cnv.info["COMMONFORWARDS"]) > 0.8
    
    def fails_meanlr2(self):
        """ checks if the MEANLR2 value is out of bounds
        """
        
        if self.cnv.genotype == "DUP":
            return float(self.cnv.info["MEANLR2"]) < 0.36
        elif self.cnv.genotype == "DEL":
            return float(self.cnv.info["MEANLR2"]) > -0.41
        
        return False
    
    def fails_no_exons(self):
        """ checks that the CNV overlaps at least one exon
        """
        
        return float(self.cnv.info["NUMBEREXONS"]) < 1
    
    def fails_frequency(self):
        """ checks that the CNV has a low population frequency.
        
        If the CNV has been observed in the unaffected controls (aka unaffected
        parents), then we can determine the population frequency, which must be
        sufficiently rare to pass.
        
        If the population frequency field is absent, assume the frequency is 0.
        """
        
        try:
            return float(self.cnv.info["ACGH_RC_FREQ50"]) > 0.01
        except KeyError:
            # If the field isn't available, assume the frequency is 0.
            return False
    
    def fails_cifer_inh(self):
        """ check that the CIFER inheritance classification isn't false_positive
        """
        
        return self.cnv.format["CIFER_INHERITANCE"] == "false_positive"
        
