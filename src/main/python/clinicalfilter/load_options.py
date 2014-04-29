""" class to load config files for clinical filtering
"""

import optparse

from clinicalfilter.load_files import *
from clinicalfilter import ped


def get_options():
    """gets the options from the command line
    """
    
    parser = optparse.OptionParser()
    parser.add_option("-p", "--ped", dest="ped_path", help="path to ped file containing cohort details for multiple trios")
    parser.add_option("-c", "--child", dest="child_path", help="path to child's VCF file")
    parser.add_option("-m", "--mother", dest="mother_path", help="path to mother's VCF file")
    parser.add_option("-f", "--father", dest="father_path", help="path to father's VCF file")
    parser.add_option("-G", "--gender", dest="child_gender", help="The child gender (male or female)")
    parser.add_option("--mom_aff_status", dest="mother_affected", help="affected status of the mother (1 = unaffected, or 2 = affected)")
    parser.add_option("--dad_aff_status", dest="father_affected", help="affected status of the father (1 = unaffected, or 2 = affected)")
    
    parser.add_option("-l", "--filter", dest="filters_path", help="path to filter file (eg filters.txt)")
    parser.add_option("-t", "--tags", dest="tags_path", help="path to tags.txt (eg tags.txt)")
    parser.add_option("--syndrome-regions", dest="cnv_regions", help="path to list of CNV regions known to occur in disorders, eg decipher_syndrome_list_20140428.txt")
    parser.add_option("--known-genes", dest="genes_path", help="path to list of known disease causative genes, eg DDG2P-reportable.txt")
    parser.add_option("--known-genes-date", dest="genes_date", help="Date that the list of known disease causative genes was last updated")
    parser.add_option("--alternate-ids", dest="alternate_ids_path", help="path to list of alternate IDs, eg personid_decipher_id_sangerid.txt")
    parser.add_option("-o", "--output", dest="output_path", default="clinical_reporting.txt", help="filename to output variant data to")
    parser.add_option("--export-vcf", dest="export_vcf", help="Folder or filename in which to export a VCF file for a proband. The script does not export a gzipped VCF file if this option is not used.")
    parser.add_option("--log", dest="loglevel", default="debug", help="level of logging to use, choose from: debug, info, warning, error or critical")
    
    (opts, args) = parser.parse_args()
    
    if opts.ped_path is None and opts.child_path is None:
        parser.error("either --ped FILENAME or --child FILENAME is required")
    if opts.ped_path and opts.child_path:
        parser.error("--ped and --child are mutually exclusive")
    if opts.filters_path is None:
        parser.error("--filter (-l) FILENAME is required")
    if opts.tags_path is None:
        parser.error("--tags (-t) FILENAME is required")
    
    return opts


class LoadOptions(object):
    """Loads the filters, weights, hierarchy and data output definitions files.
    """
    
    def set_definitions(self, opts):
        """Sets the paths to the definitions files, and parses the files for
        """
        
        self.options = opts
        
        self.filters_path = self.options.filters_path
        self.output_path = self.options.output_path
        self.tags_path = self.options.tags_path
        self.export_vcf = self.options.export_vcf
        
        self.load_definitions_files()
        self.load_trio_paths()
        
    def load_definitions_files(self):
        """loads all the config files for the script (eg filters, gene IDs)
        """
        self.filters = open_filters(self.filters_path)
        self.tags_dict = open_tags(self.tags_path)
        
        # make sure we cover all the possible ways that maximum minor allele 
        # frequencies can be named as in our VCF files. For the MAF values in 
        # the file, this should make it so that if the MAF value exists for a 
        # variant, then the variant has to pass the MAF filter
        for tag in self.tags_dict["MAX_MAF"]:
            self.filters[tag] = self.filters["MAX_MAF"]
        
        # make sure we cover all possible ways that the variants consequence ID
        # can be encoded. This should make it so that if the consequence ID 
        # exists for a variant, then that variant has to pass the consequence 
        # filter (ie have a consequence like "STOP_GAINED", 
        # "NON_SYNONYMOUS_CODING", etc)
        for tag in self.tags_dict["consequence"]:
            self.filters[tag] = self.filters["VCQ"]
        
        # if we have named a gene file, then load a dictionary of genes, and 
        # add them to the filters, so we can screen variants for being in genes 
        # known to be involved with disorders
        if self.options.genes_path is not None:
            self.known_genes = open_known_genes(self.options.genes_path)
            # include all the possible ways IDs that a gene field can be named 
            # in a VCF file
            for tag in self.tags_dict["gene"]:
                self.filters[tag] = ["list", self.known_genes]
            
            # Attempt to recover a date for the dictionary of genes
            if self.options.genes_date is not None:
                self.known_genes_date = self.options.genes_date

        else:
            self.known_genes = None
        
        # if we have named an ID mapping file, the load a dictionary of IDs and
        # alternate IDs, so we can convert between different ID schemes.
        self.ID_mapper = None
        if self.options.alternate_ids_path is not None:
            self.ID_mapper = create_person_ID_mapper(self.options.alternate_ids_path)
        
        self.cnv_regions = None
        if self.options.cnv_regions is not None:
            self.cnv_regions = open_cnv_regions(self.options.cnv_regions)
    
    def load_trio_paths(self):
        """sets the paths to the VCF files for a trio, or multiple trios.
        
        Args:
            ped_path: path to pedigree file listing all the VCF files, and their relationships
            child_path: path to VCF file for child, mutually exclusive with ped_path option
            mother_path: path to VCF file for mother, mutually exclusive with ped_path option
            father_path: path to VCF file for father, mutually exclusive with ped_path option
            childGender: gender of proband child, mutually exclusive with ped_path option
        """
        if self.options.ped_path is None:
            family = ped.Family("blank_family_ID")
            family.set_child("child", self.options.child_path, "2", self.options.child_gender)
            if self.options.mother_path is not None:
                family.set_mother("mother", self.options.mother_path, self.options.mother_affected, "2")
            if self.options.father_path is not None:
                family.set_father("father", self.options.father_path, self.options.father_affected, "1")
            
            self.families = [family]
        else:
            self.families = ped.load_families(self.options.ped_path)

