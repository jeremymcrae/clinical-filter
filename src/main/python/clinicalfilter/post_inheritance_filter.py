""" Postinheritance variant filter for clinical filtering code.

Most variants can be excluded on the basis of their individual characteristics
(minor allele frequency, quality scores etc), but some variants need to be
assessed after the standard filters
    we fail CNVs where an individual has CNVs that pass the filters on three or 
        more chromosomes (so we can only filter these out once we have found the
        CNVs that pass the filters).
    We fail CNVs with minor allele frequencies greater than 0.1%, EXCEPT for
        biallelic variants, they have a 1% threshold. We can only filter these 
        out once we have found variants that pass the different inheritance 
        models - then we can check if they are biallelelic or not.
"""

import sys
import logging

class PostInheritanceFilter(object):
    """ Apply some post inheritance filters to flagged variants
    """
    
    def __init__(self, variants):
        """intialise the class with the some definitions
        """
        
        self.variants = variants
    
    def filter_variants(self):
        """ loads trio variants, and screens for candidate variants
        """
        
        # if we have flagged CNVs on three different chroms, drop all CNVs,
        # since the sample is sufficiently anomalous
        if self.count_cnv_chroms() > 2:
            passed_vars = []
            for (var, check, inh) in self.variants:
                if not var.is_cnv():
                    passed_vars.append((var, check, inh))
            
            self.variants = passed_vars
        
        passed_vars = []
        for (var, check, inh) in self.variants:
            if inh == "Biallelelic":
                passed_vars.append((var, check, inh))
            else:
                populations = var.child.tags["MAX_MAF"]
                max_maf = var.child.find_max_allele_frequency(populations)
                print(max_maf)
                if max_maf == "NA" or float(max_maf) <= 0.1:
                    passed_vars.append((var, check, inh))
        
        self.variants = passed_vars
        
        return self.variants
    
    def count_cnv_chroms(self):
        """ count the number of different chroms that CNVs are on
        """
        
        # count the flagged CNVs
        cnv_chroms = set([])
        for (var, check, inh) in self.variants:
            if var.is_cnv():
                cnv_chroms.add(var.get_chrom())
        
        return len(cnv_chroms)

