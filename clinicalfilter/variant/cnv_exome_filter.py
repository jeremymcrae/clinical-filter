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
#        elif self.cnv.has_parents and self.fails_cifer_inh():#may be able to use self.cnv.get_has_parents and remove parents above
        elif self.fails_cifer_inh():
            passes = False
            if track_variant:
                print("failed CIFER inheritance", self.cnv.format["CIFER_INHERITANCE"])
                # elif self.fails_no_exons():
        #     passes = False
        #     if track_variant:
        #         print("failed no exons", self.info["NUMBEREXONS"])
#additional filters
#FAIL not_inherited or uncertain deletion calls if they meet at least 2 out of the 3 criteria:
#- convex_meanl2r > -1.5  
#- convex_score < 15
#- convex_mad_l2r > 0.15
#Apply total l2r on X chromosome filter of <-5000 and >7000 

        elif self.fails_x_lr2():
            passes = False
            if track_variant:
                print("fails sum mean l2r on X chromosome")

        elif self.fails_additional_filters():
            passes = False
            if track_variant:
                print("DEL fails at least 2 of 3 additional CNV filters (mean l2r < -2, score < 20, mad > 0.2)")
       
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
        if self.cnv.get_has_parents():
            return self.cnv.format["CIFER_INHERITANCE"] == "false_positive"
        else:
            return False

    def fails_x_lr2(self):
        """Fail if call is on X chromosome and sum of mean l2r on X is <-5000 
        or >7000
        """
        if self.cnv.chrom == 'X':
            x_lr2 = float(self.cnv.get_sum_x_lr2())
            if x_lr2 < -5000 or x_lr2 > 7000:
                return True
            else:
                return False
        else:
            return False

    def fails_additional_filters(self):
        """FAIL deletion calls if they meet at least 2 out of the 3 criteria:
        - convex_meanl2r > -1.5  
        - convex_score < 15
        - convex_mad_l2r > 0.15
        - filter not inherited and uncertain only
        """
        if self.cnv.genotype == "DEL":
            if self.cnv.format["CIFER_INHERITANCE"] == "not_inherited" or self.cnv.format["CIFER_INHERITANCE"] == "uncertain":
                failcount = 0
                if float(self.cnv.info["MEANLR2"]) < -1.5:
                    failcount += 1
                if float(self.cnv.info["CONVEXSCORE"]) < 15:
                    failcount += 1
                if float(self.cnv.info["MADL2R"]) > 0.15:
                    failcount += 1
                if failcount >= 2:
                    return True
                else:
                    return False
        else:#DUPs and inherited variants are OK
            return False
