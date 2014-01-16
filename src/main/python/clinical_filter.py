""" Clinical filtering for trios

Find variants in affected children that might contribute to their disorder. We
load VCF files (either named on the command line, or listed in a PED file) for
members of a family, filter for rare, functionally disruptive variants, and
assess whether each variant might affect the child's disorder. We take into 
account the parents genotypes (if available) and whether the parents are also
affected with a (the?) disorder. For variants in known disease causative genes
we check whether the inheritance patterns matches one expected for the  
inheritance models of the gene.

Usage:

python clinical_filter.py \
    --ped temp_name.ped \
    --filter filters.txt \
    --tags tags.txt \
    --known-genes known_genes.txt \
    --alternate-ids alternate_ids.txt \
    --output output_name.txt

Written by Jeremy McRae (jm33@sanger.ac.uk), derived from code by Saeed Al 
Turki (sa9@sanger.ac.uk).

"""

import sys
import optparse
import logging

from clinicalfilter import load_files
from clinicalfilter import vcf
from clinicalfilter import inheritance as inh
from clinicalfilter import ped
from clinicalfilter import reporting

class ClinicalFilter(reporting.report):
    """ filters trios for candidate variants that might contribute to a 
    probands disorder.
    """
    
    def __init__(self, defs):
        """intialise the class with the some definitions
        """
        
        # set some definitions
        self.setDefinitions(defs)
        self.first_run = True
        self.counter = 0
    
    def setDefinitions(self, defs):
        """sets the definitions for the filters and weights, as previously 
        created by the loadDefinitions class.
        """
        
        self.filters_path = defs.filters_path
        self.output_path = defs.output_path
        self.export_vcf = defs.export_vcf
        
        self.filters = defs.filters
        self.tags_dict = defs.tags_dict
        
        self.known_genes = defs.known_genes
        self.ID_mapper = defs.ID_mapper
        
        self.families = defs.families
    
    def filter_trios(self):
        """ loads trio variants, and screens for candidate variants
        """
        
        if self.first_run:
            self.printReportInputInfo(" Clinical filtering analysis")
        
        # load the trio paths into the current path setup
        for family_ID in sorted(self.families):
            self.family = self.families[family_ID]
            
            # some families have more than one child in the family, so run through each child.
            while self.family.check_all_children_analysed() == False:
                self.family.set_child()
                if self.family.child.is_affected():
                    self.vcf_loader = vcf.LoadVCFs(self.family, self.counter, len(self.families), self.filters)
                    self.variants = self.vcf_loader.get_trio_variants()
                    self.vcf_provenance = self.vcf_loader.get_vcf_provenance()
                    self.analyse_trio()
                
                self.family.set_child_examined()
            self.counter += 1
        
        # make sure we close the output file off
        if self.output_path is not None:
            self.output.close()
        
        sys.exit(0)
    
    def analyse_trio(self):
        """identify candidate variants in exome data for a single trio.
        """
        
        # report the variants that were found
        if self.first_run and self.output_path is not None:
            self.output = open(self.output_path, "w")
        
        # organise variants by gene, then find variants that fit
        # different inheritance models
        self.create_gene_dict()
        self.found_variants = []
        for gene in self.genes_dict:
            variants = self.genes_dict[gene]
            self.find_variants(variants, gene)
        
        # export the results to either tab-separated table or VCF format
        if self.output_path is not None:
            self.save_results()
        if self.export_vcf:
            self.save_vcf()
        self.first_run = False
    
    def create_gene_dict(self):
        """creates dictionary of variants indexed by gene
        """
        
        # organise the variants into entries for each gene
        self.genes_dict = {}
        for var in self.variants:
             # make sure that gene is in self.genes_dict
            if var.gene not in self.genes_dict:
                self.genes_dict[var.gene] = []
            
            # add the variant to the gene entry
            self.genes_dict[var.gene].append(var)
        
    def find_variants(self, variants, gene):
        """ finds variants that fit inheritance models
        
        Args:
            variants: list of TrioGenotype objects
            gene: gene ID as string
        """
        
        # get the inheritance for the gene (monoalleleic, biallelic, hemizygous
        # etc), but allow for times when we haven't specified a list of genes 
        # to use
        if self.known_genes is not None and gene in self.known_genes:
            gene_inheritance = self.known_genes[gene]["inheritance"]
        else:
            gene_inheritance = None
        
        # ignore intergenic variants
        if gene == None:
            return
        
        logging.debug(self.family.child.get_ID() + " " + gene + " " + str(variants) + " " + str(gene_inheritance))
        chrom_inheritance = variants[0].get_inheritance_type()
        
        if chrom_inheritance == "autosomal":
            finder = inh.Autosomal(variants, self.family, gene_inheritance)
        elif chrom_inheritance in ["XChrMale", "XChrFemale"]:
            finder = inh.Allosomal(variants, self.family, gene_inheritance)
        candidates = finder.get_candidiate_variants()
        
        for candidate in candidates:
            self.found_variants.append(candidate)


class loadDefinitions:
    """Loads the filters, weights, hierarchy and data output definitions files.
    
    We load this outside the Trios class, so that we don't have to reload the definitions each time 
    a new trio is examined.
    """
    def __init__(self, opts):
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
        """loads all the definition filters for the script (eg filters, weights, gene IDs)
        """
        self.filters = load_files.open_filters(self.filters_path)
        self.tags_dict = load_files.open_tags(self.tags_path)
        
        # make sure we cover all the possible ways that maximum minor allele frequencies can be
        # named as in our VCF files. For the MAF values in the file, this should make it so that if 
        # the MAF value exists for a variant, then the variant has to pass the MAF filter
        for tag in self.tags_dict["MAX_MAF"]:
            self.filters[tag] = self.filters["MAX_MAF"]
        
        # make sure we cover all possible ways that the variants consequence ID can be encoded. 
        # This should make it so that if the consequence ID exists for a variant, then that variant
        # has to pass the consequence filter (ie have a consequence like "STOP_GAINED", 
        # "NON_SYNONYMOUS_CODING", etc)
        for tag in self.tags_dict["consequence"]:
            self.filters[tag] = self.filters["VCQ"]
        
        # if we have named a gene file, then load a dictionary of genes, and add them to the 
        # filters, so we can screen variants for being in genes known to be involved with 
        # disorders
        if self.options.genes_path is not None:
            self.known_genes = load_files.open_known_genes(self.options.genes_path)
            # include all the possible ways IDs that a gene field can be named in a VCF file
            for tag in self.tags_dict["gene"]:
                self.filters[tag] = ["list", self.known_genes]
        else:
            self.known_genes = None
        
        # if we have named an ID mapping file, the load a dictionary of IDs and alternate IDs, so we
        # can convert between different ID schemes.
        if self.options.alternate_ids_path is not None:
            self.ID_mapper = load_files.create_person_ID_mapper(self.options.alternate_ids_path)
        else:
            self.ID_mapper = None
    
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


def get_options():
    """gets the options from the command line
    """
    
    parser = optparse.OptionParser()
    parser.add_option("-p", "--ped", dest="ped_path", help="path to ped file containing cohort details for multiple trios")
    parser.add_option("-c", "--child", dest="child_path", help="path to child's VCF file")
    parser.add_option("-m", "--mother", dest="mother_path", help="path to mother's VCF file")
    parser.add_option("-f", "--father", dest="father_path", help="path to father's VCF file")
    parser.add_option("-G", "--gender", dest="child_gender", help="The child gender (male or female)")
    parser.add_option("--mom_aff_status", dest="mother_affected", help="affected status of the mother (1 = unaffacted, or 2 = affected)")
    parser.add_option("--dad_aff_status", dest="father_affected", help="affected status of the father (1 = unaffacted, or 2 = affected)")
    
    parser.add_option("-l", "--filter", dest="filters_path", help="path to filter file (eg filters.txt)")
    parser.add_option("-t", "--tags", dest="tags_path", help="path to tags.txt (eg tags.txt)")
    parser.add_option("--known-genes", dest="genes_path", help="path to list of known disease causative genes, eg DDG2P-reportable.txt")
    parser.add_option("--alternate-ids", dest="alternate_ids_path", help="path to list of alternate IDs, eg personid_decipher_id_sangerid.txt")
    parser.add_option("-o", "--output", dest="output_path", default="clinical_reporting.txt", help="filename to output variant data to")
    parser.add_option("--export-vcf", dest="export_vcf", action="store_true", default=False, help="whether to export identified variants to a VCF file")
    parser.add_option("--log", dest="loglevel", default="debug", help="level of logging to use, choose from: debug, info, warning, error or critical")
    
    (opts, args) = parser.parse_args()
    
    if opts.ped_path and opts.child_path:
        parser.error("--ped and --child are mutually exclusive")
    if opts.filters_path is None:
        parser.error("--filter (-l) is required")
    if opts.tags_path is None:
        parser.error("--tags (-t) is required")
    
    return (opts, args)


def main():
    
    (opts, args) = get_options()
    
    # set the level of logging to generate
    loglevel = opts.loglevel
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError("Invalid log level: %s" % loglevel)
    if opts.ped_path is not None:
        log_filename = opts.ped_path + ".log"
    else:
        log_filename = "clinical-filter.log"
    logging.basicConfig(level=numeric_level, filename=log_filename)
    
    defs = loadDefinitions(opts)
    finder = ClinicalFilter(defs)
    finder.filter_trios()

if __name__ == "__main__":
    main()

