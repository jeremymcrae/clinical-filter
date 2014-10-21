""" class to load config files for clinical filtering
"""

import argparse
import sys

from clinicalfilter.load_files import open_filters, open_tags, \
    open_known_genes, create_person_ID_mapper, open_cnv_regions
from clinicalfilter import ped


def get_options():
    """gets the options from the command line
    """
    
    parser = argparse.ArgumentParser(description="Filter VCFs for inherited \
        variants in trios.")
    
    # the --ped and --child options are mutually exclusive
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ped", dest="ped", help="Path to ped file containing cohort details for multiple trios.")
    group.add_argument("--child", dest="child", help="Path to child's VCF file.")
    
    parser.add_argument("--mother", dest="mother", help="Path to mother's VCF file.")
    parser.add_argument("--father", dest="father", help="Path to father's VCF file.")
    parser.add_argument("--gender", dest="gender", help="The child's gender (male or female).")
    parser.add_argument("--mom-aff", dest="mom_aff", help="Mother's affected status (1=unaffected, or 2=affected).")
    parser.add_argument("--dad-aff", dest="dad_aff", help="Father's affected status (1=unaffected, or 2=affected).")
    
    parser.add_argument("--filter", dest="filters", required=True, help="Path to filter file (see config folder).")
    parser.add_argument( "--tags", dest="tags_path", required=True, help="Path to tags file (see config folder).")
    parser.add_argument("--syndrome-regions", dest="regions", help="Path to list of CNV regions known to occur in disorders.")
    parser.add_argument("--known-genes", dest="genes", help="Path to table of known disease causative genes.")
    parser.add_argument("--known-genes-date", dest="genes_date", help="Date that the list of known disease causative genes was last updated, used to track the version of known-genes used for analysis.")
    parser.add_argument("--alternate-ids", dest="alternate_ids", help="Path to table of alternate IDs, used to map individual IDs to their alternate study IDs.")
    parser.add_argument("-o", "--output", dest="output", help="Path for analysis output in tabular format.")
    parser.add_argument("--export-vcf", dest="export_vcf", help="Directory or file path for analysis output in VCF format.")
    parser.add_argument("--log", dest="loglevel", default="debug", help="Level of logging to use, choose from: debug, info, warning, error or critical.")
    
    # New argument added by PJ to allow DNM_PP filtering to be disabled.
    parser.add_argument("--pp-dnm-threshold", dest="pp_filter", type=float, default=0.9, help="Set PP_DNM threshold for filtering (defaults to >=0.9)")

    args = parser.parse_args()
    
    if args.child is not None and args.alternate_ids is not None:
        argparse.ArgumentParser.error("You can't specify alternate IDs when using --child")

    if args.pp_filter < 0.0 or args.pp_filter > 1:
        argparse.ArgumentParser.error("--pp-dnm-threshold must be between 0 and 1")
    
    return args


class LoadOptions(object):
    """Loads the filters, weights, hierarchy and data output definitions files.
    """
    
    def set_definitions(self, opts):
        """Sets the paths to the definitions files, and parses the files for
        """
        
        self.options = opts
        
        self.output_path = self.options.output
        self.export_vcf = self.options.export_vcf
        
        # Attempt to recover a date for the dictionary of genes
        self.known_genes_date = self.options.genes_date
        
        self.load_definitions_files()
        self.load_trio_paths()
        self.pp_filter = self.options.pp_filter
        
    def load_definitions_files(self):
        """loads all the config files for the script (eg filters, gene IDs)
        """
        self.filters = open_filters(self.options.filters)
        self.tags_dict = open_tags(self.options.tags_path)
        
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
        self.known_genes = None
        if self.options.genes is not None:
            self.known_genes = open_known_genes(self.options.genes)
            # include all the possible ways IDs that a gene field can be named 
            # in a VCF file
            for tag in self.tags_dict["gene"]:
                self.filters[tag] = ["list", self.known_genes]
        
        # if we have named an ID mapping file, the load a dictionary of IDs and
        # alternate IDs, so we can convert between different ID schemes.
        self.ID_mapper = None
        if self.options.alternate_ids is not None:
            self.ID_mapper = create_person_ID_mapper(self.options.alternate_ids)
        
        # open a list of regions associated with DECIPHER syndromes
        self.cnv_regions = None
        if self.options.regions is not None:
            self.cnv_regions = open_cnv_regions(self.options.regions)
    
    def load_trio_paths(self):
        """sets the paths to the VCF files for a trio, or multiple trios.
        """
        if self.options.ped is None:
            family = ped.Family("blank_family_ID")
            family.add_child("child", self.options.child, "2", self.options.gender)
            if self.options.mother is not None:
                family.add_mother("mother", self.options.mother, self.options.mom_aff, "2")
            if self.options.father is not None:
                family.add_father("father", self.options.father, self.options.dad_aff, "1")
            
            self.families = {family.family_id: family}
        else:
            self.families = ped.load_families(self.options.ped)

