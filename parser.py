""" class to filter VCF files, as well as sort VCF data into a gene based dictionary, and prepare
some dictionaries for exporting data. This is a bit of a mish-mash, was originally used to load and
parse VCF files. The class no longer does that, the only function left from that are the 
isPassUserFilters method, and the make_genes_dict method.
"""

import sys
import itertools
import logging


class Parser(object):
    """This class opens and parses the data from VCF files.
    """
    def __init__(self):
        pass
    
    def extractUserFields(self, raw_record):
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
        user_record = raw_record
        
        for label in self.tags_dict:
            for alternate_label in self.tags_dict[label]:
                try:
                    user_record[label] = raw_record[alternate_label]
                except KeyError:
                    continue
        
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
                print(record)
                sys.exit()
            
            # append the record to the appropriate gene entry in a dictionary
            gene_ID = record[gene_tag]
            
             # make sure that gene is in self.genes_dict
            if gene_ID not in self.genes_dict:
                self.genes_dict[gene_ID] = {}
            
            # make sure that the position is in the gene entry
            if position_key not in self.genes_dict[gene_ID]:
                self.genes_dict[gene_ID][position_key] = {}
            
            # add the VCF record to the individuals entry
            self.genes_dict[gene_ID][position_key][family_member] = record


