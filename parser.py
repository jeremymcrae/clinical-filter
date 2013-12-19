""" class to filter VCF files, as well as sort VCF data into a gene based dictionary, and prepare
some dictionaries for exporting data. This is a bit of a mish-mash, was originally used to load and
parse VCF files. The class no longer does that, the only function left from that are the 
isPassUserFilters method, and the make_genes_dict method.
"""

import sys
import itertools
import logging

import user

class Parser(object):
    """This class opens and parses the data from VCF files.
    """
    def __init__(self):
        pass
    
    def getNumber(self, values):
        '''
        This function is used for GAPI VCF files which might have multiple files
        seprated by "," at INFO columns. This should be used for generic VCF
        as the end user may not know the only the first value is returned 
        '''
        # if the string can be directly converted to a float, simply return that
        try:
            value = float(values)
        # occasionally we get comma-separated pairs (eg '.,0.639860'). Try to convert each of these 
        # in turn, if any can be converted to floats, return that value
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
        """ determines whether a value is a number, or can be converted to a number.
        
        Sometimes the MAF reported for a variant is ".", or even ".,.", which are not numbers and 
        are in fact NA values, but would cause the variant not to pass the MAF filter. instead
        check if the value can be converted to a float.
        
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
    
    def isPassUserFilters(self, record):
        """Checks whether a VCF record passes user defined criteria.
        
        Args:
            record: A dictionary entry for a single variant converted from a VCF file.
            
        Returns:
            boolean value for whether the record passes the filters
        """
        # print record
        
        show_fail_point = False
        if record['CHROM'] == '2' and record['POS'] == '38298394':
            # print record
            show_fail_point = True
        
        def show_fail():
            print key, ":", value, "not", condition, filter_values
        
        isPass = True
        for key, value in record.items():
            # # don't bother to check the value if we haven't got one recorded for the key
            # if value is None:
            #     continue
            if key in self.filters:
                condition = self.filters[key][0]
                filter_values = self.filters[key][1]
                if condition == "list":
                    if value not in filter_values:
                        isPass = False
                        if show_fail_point:
                            show_fail()
                        break   
                elif condition == "greater_than":
                    value = self.getNumber(value)
                    if value < filter_values and self.is_number(value):
                        isPass = False
                        if show_fail_point:
                            show_fail()
                        break
                elif condition == "smaller_than":
                    # some of the MAF values are 1 - MAF due to being for a population that was 
                    # genotyped on the opposing strand. We need to convert those back.
                    # if key in self.tags_dict["MAX_MAF"]:
                    #     value = self.getNumber(value)
                    #     if self.is_number(value):
                    #         if value > 0.5:
                    #             value = 1 - value
                    value = self.getNumber(value)
                    if value > filter_values and self.is_number(value):
                        isPass = False
                        if show_fail_point:
                            show_fail()
                        break
                elif condition == "equal":
                    if value != filter_values:
                        isPass = False
                        if show_fail_point:
                            show_fail()
                        break
                elif condition == "not":
                    if value == filter_values:
                        isPass = False
                        if show_fail_point:
                            show_fail()
                        break
                elif condition == "multiple_not":
                    # sometimes we want to make sure two related keys don't contain certain values.
                    # Specifically, we want to exclude variants where polyphen annotation is 
                    # "benign" and sift annotation is "tolerated". These are provided as a list of 
                    # tuples. We catch the first key, and then check the second key at the same time
                    has_all_values = True
                    for key, filter_value in filter_values:
                        # pull out the value for each key in the list, and check whether it contains
                        # the filter value
                        temp_value = record[key]
                        if temp_value is None:
                            has_all_values = False
                        # the value should be something like benign(0.05), or tolerated(0.3), and we
                        # just check if "benign", or "tolerated" are in the corresponding values
                        elif filter_value not in temp_value:
                            has_all_values = False
                    if has_all_values:
                        isPass = False
                        if show_fail_point:
                            show_fail()
                        break
                elif condition in ["startswith", "endswith"]:
                    if condition == "startswith":
                        if not value.startswith(filter_values):
                            isPass = False
                            break
                    elif condition == "endswith":
                        if not value.endswith(filter_values):
                            isPass = False
                            break
                elif condition == "range":
                    start, end = filter_values
                    if self.getNumber(value) < start or self.getNumber(value) > end:
                        isPass = False
                        break
        return isPass
    
    def find_max_allele_frequency(self, record):
        """gets the maximum allele frequency for a variant in a VCF record.
        
        Finds the maximum allele frequency recorded for a variant iacross different populations.
        
        Args:
            record: A VCF dictionary record for a single variant.
            
        Returns:
            the maximum allele frequency found within the populations recorded in the variant record
        """
        
        max_allele_frequency = -100
        # run through all the possible populations in the VCF file (typically the 1000 Genomes 
        # populations (AFR_AF, EUR_AF etc), an internal popuation (DDD_AF), and a AF_MAX field)
        for key in self.tags_dict["MAX_MAF"]:
            if key in record:
                number = self.getNumber(record[key])
                if not self.is_number(number):
                    continue
                # if number > 0.5:
                #     number = 1 - number
                #     record[key] = str(number)
                if number > max_allele_frequency:
                    max_allele_frequency = number
        
        # return NA for variants without MAF recorded
        if max_allele_frequency == -100:
            max_allele_frequency = "NA"
        
        return str(max_allele_frequency)
    
    def extractUserFields(self, raw_record, raw_header):
        '''Selects the data to report as output.
        
        If the user requested certain fields for the output , then extract them. If not, print 
        everything but in the same order of the original VCF record.
        
        Called by var_under_models()
        
        Args:
            raw_record: the VCF record for a variant for an individuals
            raw_header: the parsed header position of the VCF file for the individuals
        
        Returns:
            user_record: 
        '''
        user_record = {}
        if self.orders['output'] != "" and self.orders['output'] != None or self.orders['output'] == "*":
            for label in self.orders["output"]: 
                # allow for VCF files with different IDs, and make sure we still catch the right
                # data, by checking the different tags in self.tags_dict
                if label in self.tags_dict:
                    for alternate_label in self.tags_dict[label]:
                        try:
                            user_record[label] = raw_record[alternate_label]
                        except KeyError:
                            continue
                else:  
                    user_record[label] = raw_record[label]
        else:
            for label in raw_header:
                user_record[label] = raw_record[label]
        return user_record
    
    def make_genes_dict(self, vcf, family_member):
        """Parses variant data into gene based dictionary entries.
        
        Args:
            vcf: Dictionary entry for individual from VCF file
            family_member: individuals position in family trio eg 'mother', 'father', 'child'
            
        Returns:
            Nothing, instead creates a dictionary belonging to the parser class.
            
        """
        for key in vcf['data']:
            position_key = key[1] # extract pos
            record = vcf['data'][key]
            
            # pull out a gene name for the variant, no matter how the gene field is named
            gene_tag = None
            for tag in self.tags_dict["gene"]:
                if tag in record:
                    gene_tag = tag
            
            if gene_tag is None:
                print record
                sys.exit()
            
            # append the record to the appropriate gene entry in a dictionary
            gene_ID = record[gene_tag]
            
             # make sure that gene is in self.genes_dict
            if gene_ID not in self.genes_dict:
                self.genes_dict[gene_ID] = {"positions": {}}
            
            # make sure that the position is in the gene entry
            if position_key not in self.genes_dict[gene_ID]["positions"]:
                self.genes_dict[gene_ID]["positions"][position_key] = {}
            
            # make sure that the individual is in the variant entry
            if family_member not in self.genes_dict[gene_ID]["positions"][position_key]:
                self.genes_dict[gene_ID]["positions"][position_key][family_member] = {}
            
            # add the VCF record to the individuals entry
            self.genes_dict[gene_ID]["positions"][position_key][family_member] = record
    
    # def getWeights(self, record):
    #     """Gets the weights recorded for a single variant entry for an individual.
        
    #     Args:
    #         record: a VCF dictionary record for an individual for a single variant
        
    #     Returns:
    #         results: a list of weight values, sorted by the same order as in self.orders['weights']
    #     """

    #     weights_dict = {}
    #     # fill the weights_dict by default with zero weights, this covers all the possible weight
    #     # keys, even if they don't exist in the VCF variant record that is being processed
    #     for key in self.weights:
    #         weights_dict[key] = 0.0
        
    #     for tag_name in self.weights:
    #         # ignore weights that aren't listed in the VCF record
    #         if tag_name not in record:
    #             continue
            
    #         tag_value = record[tag_name]
    #         for condition in self.weights[tag_name]:
    #             if condition == "list":
    #                 for lst in self.weights[tag_name]["list"]:
    #                     if tag_value in lst[0]:
    #                         weights_dict[tag_name] = float(lst[1])
    #             elif condition in ["greater_than","smaller_than", "equal"]:
    #                 try:
    #                     tag_value = float(tag_value)
    #                     cutoff = self.weights[tag_name][condition][0][0]
    #                     wt = self.weights[tag_name][condition][0][1]
    #                 except:
    #                     print tag_name, condition, self.weights[tag_name][condition]
    #                     print "Either filter cutoff %s of VCF file %s is not numerical."  % (cutoff,tag_value)
    #                     print "Please correct your filters or / and VCF file. The program will stop now."
    #                     sys.exit(0)
    #                 if condition == "greater_than":
    #                     if tag_value > cutoff:
    #                         weights_dict[tag_name] = wt
    #                 elif condition == "smaller_than":
    #                     if tag_value < cutoff:
    #                         weights_dict[tag_name] = wt
    #                 elif condition == "equal":
    #                     if tag_value == cutoff:
    #                         weights_dict[tag_name] = wt
    #             elif condition in ["startswith","endswith"]:
    #                 cutoff = self.weights[tag_name][condition][0][0]
    #                 wt = self.weights[tag_name][condition][0][1]
    #                 if condition == "startswith":
    #                     if tag_value.startswith(cutoff):
    #                         weights_dict[tag_name] = wt
    #                 elif condition == "endswith":
    #                     if tag_value.endswith(cutoff):
    #                             weights_dict[tag_name] = wt
    #             elif condition == "range":
    #                 start, end = self.weights[tag_name][condition][0][0]
    #                 wt = self.weights[tag_name][condition][0][1]
    #                 tag_value = float(tag_value)
    #                 if tag_value  > start and tag_value < end:
    #                     weights_dict[tag_name] = wt
        
    #     # and standardise the weight values to the list found in self.orders['weights']
    #     weights = []
    #     for k in self.orders["weights"]:
    #         weights.append(weights_dict[k])
    #     return weights
    
    def construct_holder(self, type="numerical"):
        '''
        Construct a holder dictionary where keys are either numeric or strings.
        Numerical keys will be used for filtering.
        String keys will be used to print results and for the sunburst chart labels.
        
        I don't really understand what this is doing, or why. Currently it's used in two places, in 
        get_var_under_models(), and later in getLabels(), which is called by printStatsTable(), 
        printResults() and save2R(), but the results aren't actually used by those functions, so I 
        have commented out those lines. 
        '''
        
        # do we want to use names or not
        type_dict = {"numerical": False, "labels": True}
        names = type_dict[type]
        mydict = self.get_tree_lvls(names)
        
        # I believe this extracts all the levels from the weights filtering criteria
        lvls_lst = []
        for lvl in sorted(mydict.keys()):
            lvls_lst.append(sorted(mydict[lvl]))
        
        keys =  tuple(itertools.product(*lvls_lst)) # make permutation of all possible keys
        holder = {}
        for key in keys:
            if key not in holder:
                holder[key] = {}
        
        return holder
    
    def get_tree_lvls(self, names=False):
        ''' Creates a dictionary of entries from the weight filter criteria.
        
        Example self.AR_results = {LOF:GERP:CONDEL:"empty_space_for_raw_VCF_line"}
        where LOF, GERP, CONDEL are int (0,1,2,..,n)
        
        NOTE: I don't really understand this function, the how and why.
        '''
        
        levels = {}
        # get the weight from each weight class and all its conditions
        for i, wt_class in enumerate(self.orders["weights"]):
            for condition in self.weights[wt_class]:
                for lst in self.weights[wt_class][condition]:
                    
                    # change the values depending on whether we include names or not
                    if names == True:
                        name = lst[2] # the weight value
                        weight = wt_class
                    else:
                        name = lst[1]
                        weight = i
                    
                    if weight not in levels:
                        levels[weight] = []
                    levels[weight].append(name)
                        
        return levels



