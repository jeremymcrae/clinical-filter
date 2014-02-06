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

from clinicalfilter import vcf
from clinicalfilter.inheritance import *
from clinicalfilter.reporting import Report
from clinicalfilter.load_options import LoadOptions, get_options

class ClinicalFilter(LoadOptions, Report):
    """ filters trios for candidate variants that might contribute to a 
    probands disorder.
    """
    
    def __init__(self, opts):
        """intialise the class with the some definitions
        """
        
        # set some definitions
        self.set_definitions(opts)
        self.first_run = True
        self.counter = 0
    
    def filter_trios(self):
        """ loads trio variants, and screens for candidate variants
        """
        
        if self.first_run:
            self.printReportInputInfo(" Clinical filtering analysis")
        
        # load the trio paths into the current path setup
        for family_ID in sorted(self.families):
            self.family = self.families[family_ID]
            
            # some families have more than one child in the family, so run 
            # through each child.
            self.family.set_child()
            while self.family.child is not None:
                if self.family.child.is_affected():
                    self.vcf_loader = vcf.LoadVCFs(self.family, self.counter, \
                        len(self.families), self.filters, self.tags_dict)
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
            self.found_variants += self.find_variants(variants, gene)
        
        # export the results to either tab-separated table or VCF format
        if self.output_path is not None:
            self.save_results()
        if self.export_vcf is not None:
            self.save_vcf()
        self.first_run = False
    
    def create_gene_dict(self):
        """creates dictionary of variants indexed by gene
        """
        
        # organise the variants into entries for each gene
        self.genes_dict = {}
        for var in self.variants:
             # make sure that gene is in self.genes_dict
            if var.get_gene() not in self.genes_dict:
                self.genes_dict[var.get_gene()] = []
            
            # add the variant to the gene entry
            self.genes_dict[var.get_gene()].append(var)
        
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
            gene_inh = self.known_genes[gene]["inheritance"]
        else:
            gene_inh = None
        
        # ignore intergenic variants
        if gene is None:
            return []
        
        logging.debug(self.family.child.get_ID() + " " + gene + " " + \
            str(variants) + " " + str(gene_inh))
        chrom_inheritance = variants[0].get_inheritance_type()
        
        if chrom_inheritance == "autosomal":
            finder = Autosomal(variants, self.family, self.known_genes, gene_inh)
        elif chrom_inheritance in ["XChrMale", "XChrFemale", "YChrMale"]:
            finder = Allosomal(variants, self.family, self.known_genes, gene_inh)
        
        return finder.get_candidate_variants()


def main():
    
    options = get_options()
    
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
    
    finder = ClinicalFilter(options)
    finder.filter_trios()

if __name__ == "__main__":
    main()

