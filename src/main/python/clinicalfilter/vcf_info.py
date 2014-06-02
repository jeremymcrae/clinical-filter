""" class for filtering SNVs based on VCF INFO fields
"""

import copy

class VcfInfo(object):
    """ parses the VCF info field, and checks whether the variant passes 
    filtering criteria.
    """
    
    show_fail_point = False
    
    def add_info(self, info_values, tags):
        """Parses the INFO column from VCF files.
        
        Args:
            info_values: INFO text from a line in a VCF file
            tags: the tags dict
        """
        
        self.tags = tags
        self.info = {}
        
        for item in info_values.split(";"):
            if "=" in item:
                try:
                    key, value = item.split("=")
                except ValueError:
                    pos = item.index("=")
                    key = item[:pos]
                    value = item[pos + 1:]
            else:
                key, value = item, True
            self.info[key] = value
        
        # add the filter value, as we filter with the info dict
        self.info["FILTER"] = self.filter
        
        self.add_gene_from_info()
        self.add_consequence()
    
    def has_info(self):
        """ checks if the INFO field has been parsed and added to the object
        """
        
        return hasattr(self, "info")
    
    def add_gene_from_info(self):
        """ adds a gene to the var using the info. CNVs and SNVs act differently
        """
        
        # sometimes the variant lacks an HGNC field
        if "HGNC" not in self.info:
            if "gene" not in self.info:
                self.gene = None
            else:
                self.gene = self.info["gene"]
        else:
            self.gene = self.info["HGNC"]
    
    def add_consequence(self):
        """ makes sure a consequence field is available in the info dict
        """
        
        # for consequence in self.tags["consequence"]:
        #     if consequence in self.info:
        #         self.info["CQ"] = self.info[consequence]
        
        if "CQ" not in self.info:
            self.info["CQ"] = None
            
    def is_lof(self):
        """ checks if a variant has a loss-of-function consequence
        """
        
        # define the set of loss-of-function consequences
        lof_consequences = set(["transcript_ablation","splice_donor_variant", \
            "splice_acceptor_variant", "frameshift_variant", "stop_gained", \
            "coding_sequence_variant"])
        
        return self.info["CQ"] in lof_consequences
    
    def get_number(self, values):
        """ converts a string into a number
        
        This function is used for GAPI VCF files which might have multiple files
        seprated by "," at INFO columns. This should be used for generic VCF
        as the end user may not know the only the first value is returned 
        """
        # if the string can be directly converted to a float, simply return that
        try:
            value = float(values)
        # occasionally we get comma-separated pairs (eg '.,0.639860'). Try to 
        # convert each of these in turn, if any can be converted to floats, 
        # return that value
        except ValueError:
            values = values.split(",")
            for value in values:
                try:
                    value = float(value)
                    break
                except ValueError:
                    pass
        except:
            value = values
        return value
    
    def is_number(self, value):
        """ determines whether a value represents a number.
        
        Sometimes the MAF reported for a variant is ".", or even ".,.", which 
        are not numbers and are in fact NA values, but would cause the variant
        not to pass the MAF filter. instead check if the value can be 
        converted to a float.
        
        Args:
            value: a string or other number
        
        Returns:
            True or False for whether the value can be converted to a float.
        """
        
        if value is None:
            return False
        
        try:
            value = float(value)
            return True
        except ValueError:
            return False
        
        return False
    
    def find_max_allele_frequency(self, populations):
        """gets the maximum allele frequency for a variant in a VCF record
        
        Finds the maximum allele frequency recorded for a variant across
        different populations.
        
        Args:
            populations: list of population IDs to search
          
        Returns:
            the maximum allele frequency found within the populations in the
            variant record
        """
        
        max_allele_frequency = -100
        # run through all the possible populations in the VCF file (typically 
        # the 1000 Genomes populations (AFR_AF, EUR_AF etc), an internal 
        # popuation (DDD_AF), and a AF_MAX field)
        for key in populations:
            if key in self.info:
                number = self.get_number(self.info[key])
                if not self.is_number(number):
                    continue
                # if number > 0.5:
                #     number = 1 - number
                if number > max_allele_frequency:
                    max_allele_frequency = number
        
        # return NA for variants without MAF recorded
        if max_allele_frequency == -100:
            max_allele_frequency = "NA"
        
        return str(max_allele_frequency)
    
    def show_fail(self, key, value, condition, filter_values):
        """ prints why a named variant has failed filtering
        """
        print(str(key) + ": " + str(value) + " not " + str(condition) + \
                  " " + str(filter_values))
    
    def passes_filters(self, filters):
        """Checks whether a VCF record passes user defined criteria.
        
        Args:
            filters: A dictionary of filtering criteria.
            
        Returns:
            boolean value for whether the variant passes the filters
        """
        
        self.show_fail_point = False
        if self.get_chrom() == "X" and self.get_position() == "48935752":
            self.show_fail_point = True
        
        passes = True
        for key in self.info:
            if key not in filters:
                continue
            
            value = self.info[key]
            condition = filters[key][0]
            filter_values = filters[key][1]
            
            if condition == "list":
                passes = self.passes_list(value, filter_values)
            elif condition == "smaller_than":
                passes = self.passes_smaller_than(value, filter_values)
                
            if passes == False:
                break
        
        if passes == False and self.show_fail_point:
            self.show_fail(key, value, condition, filter_values)
        
        return passes
    
    def passes_list(self, value, filter_values):
        """ checks whether the vcf value is within a list 
        """
        
        if filter_values is None:
            return False
        
        return value in filter_values
    
    def passes_smaller_than(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        # some of the MAF values are 1 - MAF due to being for a population that 
        # was genotyped on the opposing strand. We need to convert those back.
        # if key in self.tags["MAX_MAF"]:
        #     value = self.get_number(value)
        #     if self.is_number(value):
        #         if value > 0.5:
        #             value = 1 - value
        value = self.get_number(value)
        try:
            value > filter_values
        except TypeError:
            return True
        
        return value <= filter_values and self.is_number(value)
