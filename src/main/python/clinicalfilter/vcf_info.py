""" class for filtering SNVs based on VCF INFO fields
"""

class VcfInfo(object):
    """ parses the VCF info field, and checks whether the variant passes 
    filtering criteria.
    """
    
    # Here are the VEP consequences, ranked in severity from the most severe to
    # the least severe as defined at: 
    # http://www.ensembl.org/info/genome/variation/predicted_data.html
    severity = {"transcript_ablation": 0, "splice_donor_variant": 1, \
        "splice_acceptor_variant": 2, "stop_gained": 3, "frameshift_variant": 4, \
        "stop_lost": 5, "initiator_codon_variant": 6, "inframe_insertion": 7, \
        "inframe_deletion": 8, "missense_variant": 9, \
        "transcript_amplification": 10, "splice_region_variant": 11, \
        "incomplete_terminal_codon_variant": 12, "synonymous_variant": 13, \
        "stop_retained_variant": 14, "coding_sequence_variant": 15, \
        "mature_miRNA_variant": 16, "5_prime_UTR_variant": 17, \
        "3_prime_UTR_variant": 18, "intron_variant": 19, \
        "NMD_transcript_variant": 20, "non_coding_exon_variant": 21, \
        "nc_transcript_variant": 22, "upstream_gene_variant": 23, \
        "downstream_gene_variant": 24, "TFBS_ablation": 25, \
        "TFBS_amplification": 26, "TF_binding_site_variant": 27, \
        "regulatory_region_variant": 28, "regulatory_region_ablation": 29, \
        "regulatory_region_amplification": 30, "feature_elongation": 31, \
        "feature_truncation": 32, "intergenic_variant": 33}
    
    debug_chrom = None
    debug_pos = None
    
    def add_info(self, info_values, tags):
        """Parses the INFO column from VCF files.
        
        Args:
            info_values: INFO text from a line in a VCF file
            tags: the tags dict
        """
        
        self.tags = tags
        
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
        
        self.set_consequence()
        self.add_gene_from_info()
    
    def has_info(self):
        """ checks if the INFO field has been parsed and added to the object
        """
        
        return self.info != {}
    
    def add_gene_from_info(self):
        """ adds a gene to the var using the info. CNVs and SNVs act differently
        """
        
        # sometimes the variant lacks an HGNC field
        if "HGNC" in self.info:
            self.gene = self.info["HGNC"]
        else:
            self.gene = None
            for gene_tag in self.tags["gene"]:
                if gene_tag in self.info:
                    self.gene = self.info[gene_tag]
                    break
    
    def set_consequence(self):
        """ makes sure a consequence field is available in the info dict
        """
        
        if "CQ" not in self.info:
            self.info["CQ"] = None
        
        # some variants have multiple alts, so we need to select the alt with 
        # the most severe consequence. However, in at least one version of the 
        # VCFs, one of the alts could have zero depth, which I believe resulted 
        # from the population based multi-sample calling. We need to drop the
        # consequences recorded for zero-depth alternate alleles before finding
        # the most severe.
        if "," in self.alt_allele:
            # drop the consequence for zero-depth alleles
            if "AC" in self.info:
                cq = self.info["CQ"].split(",")
                
                # check if alleles are present in the individual
                counts = [ x != "0" for x in self.info["AC"].split(",") ]
                
                # exclude the consequences for alleles with zero alleles
                cq[:] = [ item for i,item in enumerate(cq) if counts[i] ]
                
                # convert the alleles back into a comma separated list
                self.info["CQ"] = ",".join(cq)
                
                if "HGNC" in self.info:
                    enst = self.info["ENST"].split(",")
                    hgnc = self.info["HGNC"].split(",")
                    hgnc[:] = [ item for i,item in enumerate(hgnc) if counts[i] ]
                    enst[:] = [ item for i,item in enumerate(enst) if counts[i] ]
                    self.info["HGNC"] = ",".join(hgnc)
                    self.info["ENST"] = ",".join(enst)
            
            consequences = self.info["CQ"].split(",")
            self.info["CQ"] = self.get_most_severe_consequence(consequences)
            
            if "HGNC" in self.info:
                # and get the unique HGNC and transcripts
                self.info["HGNC"] = ",".join(sorted(set(self.info["HGNC"].split(","))))
                self.info["ENST"] = ",".join(sorted(set(self.info["ENST"].split(","))))
             
    def get_most_severe_consequence(self, consequences):
        """ get the most severe consequence from a list of vep consequence terms
        
        Args:
            consequences: list of VEP consequence strings
        
        Returns:
            the most severe consequence string
        """
        
        most_severe = ""
        most_severe_score = 1000
        for cq in self.severity:
            if self.severity[cq] < most_severe_score:
                most_severe = cq
                most_severe_score = self.severity[cq]
        
        return most_severe
       
    def is_lof(self):
        """ checks if a variant has a loss-of-function consequence
        """
        
        # define the set of loss-of-function consequences
        lof_consequences = set(["transcript_ablation", "splice_donor_variant", \
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
    
    def passes_filters(self, filters):
        """Checks whether a VCF record passes user defined criteria.
        
        Args:
            filters: A dictionary of filtering criteria.
            
        Returns:
            boolean value for whether the variant passes the filters
        """
        
        pass_value, key = self.check_filters(filters)
        
        return pass_value
    
    def check_filters(self, filters, debug=False):
        """Checks whether a VCF record passes user defined criteria.
        
        Args:
            filters: A dictionary of filtering criteria.
            debug: if True, stop at the first failed filter, if false, 
                check all the filters.
            
        Returns:
            boolean value for whether the variant passes the filters
        """
        
        pass_value = True
        for key in self.info:
            if key not in filters:
                continue
            
            value = self.info[key]
            condition = filters[key][0]
            filter_values = filters[key][1]
            
            if condition == "list":
                pass_value = self.passes_list(value, filter_values)
            elif condition == "smaller_than":
                pass_value = self.passes_smaller_than(value, filter_values)
            
            # stop at the first failed filter
            if not pass_value:
                break
                
        return pass_value, key
    
    def passes_filters_with_debug(self, filters):
        """Checks whether a VCF record passes user defined criteria.
        
        This method replaces passes_filters() when we specify a chromosome and
        position for debugging the filtering.
        
        Args:
            filters: A dictionary of filtering criteria.
            
        Returns:
            boolean value for whether the variant passes the filters
        """
        
        pass_value, key = self.check_filters(filters)
        
        if pass_value == False and self.get_position() == self.debug_pos:
            value = self.info[key]
            condition = filters[key][0]
            filter_values = filters[key][1]
            
            print("failed {0}: {1} not {2} {3}".format(key, value, condition, \
                filter_values))
        
        return pass_value
    
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
        
        if not self.is_number(value):
            return True
        
        return value <= filter_values
