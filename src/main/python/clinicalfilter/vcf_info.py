""" class for filtering SNVs based on VCF INFO fields
"""


class VcfInfo(object):
    """ parses the VCF info field, and checks whether the variant passes 
    filtering criteria.
    """
    
    def __init__(self):
        """
        """
        
        self.show_fail_point = False
    
    def add_info(self, info_values):
        """Parses the INFO column from VCF files.
        
        Args:
            info_values: INFO text from a line in a VCF file
        """
        
        self.info = {}
        
        for item in info_values.split(";"):
            if "=" in item:
                key, value = item.split("=")
            else:
                key, value = item, True
            self.info[key] = value
        
        # add the filter value, as we filter with the info dict
        self.info["FILTER"] = self.filter
        
        # sometimes the variant lacks an HGNC field, but does have a HGNC_ALL
        # entry. 
        if "HGNC" not in self.info and "HGNC_ALL" in self.info:
            # if there is only one gene for HGNC_ALL, just use that as the gene
            self.info["HGNC"] = self.info["HGNC_ALL"]
            # if len(self.info["HGNC_ALL"].split(",")) == 1:
            # # if there is more than one, call the gene NA, as CNVs can encompass
            # # multiple genes, but we don't want those to match anything in the 
            # # DDG2P list.
            # else:
            #     self.info["HGNC"] = "NA"
        elif "HGNC" not in self.info and "HGNC_ALL"  not in self.info and "NUMBERGENES" in self.info:
            self.info["HGNC"] = None
            if int(self.info["NUMBERGENES"]) > 0:
                self.info["HGNC"] = "."
        # make sure we have gene and consequence keys in the info dict, for 
        # the filter to work with
        elif "HGNC" not in self.info:
            self.info["HGNC"] = None
        
        self.gene = self.info["HGNC"]
        
        if "CQ" not in self.info:
            self.info["CQ"] = None
    
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
    
    def show_fail(self, key, value, condition, filter_values):
        print(str(key) + ": " + str(value) + " not " + str(condition) + \
                  " " + str(filter_values))
    
    def passes_filters(self, filters):
        """Checks whether a VCF record passes user defined criteria.
        
        Args:
            record: A dictionary entry for a single variant converted from a
            VCF file.
            
        Returns:
            boolean value for whether the record passes the filters
        """
        
        self.show_fail_point = False
        if self.chrom == '1' and self.position == '18937':
            show_fail_point = True
        
        passes = True
        for key in self.info:
            if key not in filters:
                continue
            
            value = self.info[key]
            condition = filters[key][0]
            filter_values = filters[key][1]
            
            if condition == "list" and self.fails_list(value, filter_values):
                passes = False
                break
            elif condition == "greater_than" and self.fails_greater_than(value, filter_values):
                passes = False
                break
            elif condition == "smaller_than" and self.fails_smaller_than(value, filter_values):
                passes = False
                break
            elif condition == "equal" and self.fails_equal(value, filter_values):
                passes = False
                break
            elif condition == "not" and self.fails_not(value, filter_values):
                passes = False
                break
            elif condition == "multiple_not" and self.fails_multiple_not(value, filter_values):
                passes = False
                break
            elif condition == "startswith" and self.fails_start_string(value, filter_values):
                passes = False
                break
            elif condition == "endswith" and self.fails_end_string(value, filter_values):
                passes = False
                break
            elif condition == "range" and self.fails_range(value, filter_values):
                passes = False
                break
                
            if passes == False:
                break
        
        if passes == False:
            if self.show_fail_point:
                self.show_fail(key, value, condition, filter_values)
        
        return passes
    
    def fails_list(self, value, filter_values):
        """ checks whether the vcf value is within a list 
        """
        
        fails = False
        if value not in filter_values:
            fails = True
        
        return fails
    
    def fails_greater_than(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        value = self.get_number(value)
        try:
            value > filter_values
        except TypeError:
            return False
            
        fails = False
        if value < filter_values and self.is_number(value):
            fails = True
        
        return fails
    
    def fails_smaller_than(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        # some of the MAF values are 1 - MAF due to being for a population that 
        # was genotyped on the opposing strand. We need to convert those back.
        # if key in self.tags_dict["MAX_MAF"]:
        #     value = self.get_number(value)
        #     if self.is_number(value):
        #         if value > 0.5:
        #             value = 1 - value
        value = self.get_number(value)
        try:
            value > filter_values
        except TypeError:
            return False
        
        fails = False
        if value > filter_values and self.is_number(value):
            fails = True
        
        return fails
                    
    def fails_equal(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        fails = False
        if value != filter_values:
            fails = True
            
        return fails
    
    def fails_not(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        fails = False
        if value == filter_values:
            fails = True
            
        return fails
    
    def fails_range(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        fails = False
        start, end = filter_values
        if self.get_number(value) < start or self.get_number(value) > end:
            fails = True
            
        return fails
    
    def fails_start_string(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        fails = False
        if not value.startswith(filter_values):
            fails = True
        
        return fails
        
    def fails_end_string(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        fails = False
        if not value.endswith(filter_values):
            fails = True
        
        return fails
        
    def fails_multiple_not(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        # sometimes we want to make sure two related keys don't contain certain
        # values. Specifically, we want to exclude variants where polyphen 
        # annotation is "benign" and sift annotation is "tolerated". These are
        # provided as a list of  tuples. We catch the first key, and then check
        # the second key at the same time
        has_all_values = True
        for key, filter_value in filter_values:
            # pull out the value for each key in the list, and check whether it
            # contains the filter value
            temp_value = self.info[key]
            if temp_value is None:
                has_all_values = False
            # the value should be something like benign(0.05), or 
            # tolerated(0.3),and we just check if "benign", or "tolerated" are 
            # in the corresponding values
            elif filter_value not in temp_value:
                has_all_values = False
        
        fails = False
        if has_all_values:
            passes = True
        
        return fails

