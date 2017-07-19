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
        elif self.fails_meanlr2():
            passes = False
            if track_variant:
                print("failed meanlr2", self.cnv.info["MEANLR2"])
        elif self.fails_commmon_forwards():
            passes = False
            if track_variant:
                print("failed commonforwards", self.cnv.info["COMMONFORWARDS"])
        elif self.fails_cifer_inh():
            passes = False
            if track_variant:
                print("failed CIFER inheritance", self.cnv.format["CIFER_INHERITANCE"])
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
            return abs(float(self.cnv.info["MEANLR2"])/float(self.cnv.info["MADL2R"])) < 10
        except ZeroDivisionError:
            return True
    
    def fails_meanlr2(self):
        """ checks if the MEANLR2 value is out of bounds
        """
        
        if self.cnv.genotype == "DUP":
            return float(self.cnv.info["MEANLR2"]) < 0.4
        elif self.cnv.genotype == "DEL":
            return float(self.cnv.info["MEANLR2"]) > -0.5
        
        return False
    
    def fails_commmon_forwards(self):
        """ checks if the COMMONFORWARDS value is too high
        """
        
        return float(self.cnv.info["COMMONFORWARDS"]) > 0.8
    
    def fails_no_exons(self):
        """ checks that the CNV overlaps at least one exon
        """
        
        return float(self.cnv.info["NUMBEREXONS"]) < 1
    
    def fails_cifer_inh(self):
        """ check that the CIFER inheritance classification isn't false_positive
        """
        
        return self.cnv.format["CIFER_INHERITANCE"] == "false_positive"
