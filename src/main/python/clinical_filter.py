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
    --output output_name.txt \
    --pp-dnm-threshold threshold_as_float (default 0.9)

Written by Jeremy McRae (jm33@sanger.ac.uk), derived from code by Saeed Al
Turki (sa9@sanger.ac.uk) and Jeff Barrett.
"""

import sys
import logging

from clinicalfilter.load_vcfs import LoadVCFs
from clinicalfilter.inheritance import Allosomal, Autosomal, CNVInheritance
from clinicalfilter.post_inheritance_filter import PostInheritanceFilter
from clinicalfilter.reporting import Report
from clinicalfilter.load_options import LoadOptions, get_options

class ClinicalFilter(LoadOptions):
    """ filters trios for candidate variants that might contribute to a
    probands disorder.
    """
    
    def __init__(self, opts):
        """intialise the class with the some definitions
        """
        
        self.set_definitions(opts)
        self.report = Report(self.output_path, self.export_vcf, self.ID_mapper,
            self.known_genes_date)
    
    def filter_trios(self):
        """ loads trio variants, and screens for candidate variants
        """
        
        self.vcf_loader = LoadVCFs(len(self.families), self.known_genes, \
            self.debug_chrom, self.debug_pos, self.last_base)
        
        # load the trio paths into the current path setup
        for family_ID in sorted(self.families):
            self.family = self.families[family_ID]
            
            # some families have more than one child in the family, so run
            # through each child.
            self.family.set_child()
            while self.family.child is not None:
                if self.family.child.is_affected():
                    variants = self.vcf_loader.get_trio_variants(self.family, self.pp_filter)
                    self.vcf_provenance = self.vcf_loader.get_trio_provenance()
                    self.analyse_trio(variants)
                
                self.family.set_child_examined()
        
        sys.exit(0)
    
    def analyse_trio(self, variants):
        """identify candidate variants in exome data for a single trio.
        
        takes variants that passed the initial filtering from VCF loading, and
        splits the variants into groups for each gene with variants. Then
        analyses variants in a single gene (so we can utilise the appropriate
        inheritance mechanisms for that gene), before running some
        pos-inheritance filters, and exporting the data (ir required).
        
        Args:
            variants: list of TrioGenotypes objects
        """
        
        # organise variants by gene, then find variants that fit
        # different inheritance models
        genes_dict = self.create_gene_dict(variants)
        found_vars = []
        for gene in genes_dict:
            gene_vars = genes_dict[gene]
            found_vars += self.find_variants(gene_vars, gene)
        
        # remove any duplicate variants (which might ocur due to CNVs being
        # checked against all the genes that they encompass)
        found_vars = self.exclude_duplicates(found_vars)
        
        # apply some final filters to the flagged variants
        post_filter = PostInheritanceFilter(found_vars, self.family, self.debug_chrom, self.debug_pos)
        found_vars = post_filter.filter_variants()
        
        # export the results to either tab-separated table or VCF format
        self.report.export_data(found_vars, self.family, \
            self.vcf_loader.child_header, self.vcf_provenance)
    
    def create_gene_dict(self, variants):
        """creates dictionary of variants indexed by gene
        
        Args:
            variants: list of TrioGenotypes objects
        
        Returns:
            dictionary of variants indexed by HGNC symbols
        """
        
        # organise the variants into entries for each gene
        genes = {}
        for var in variants:
            # variants (particularly CNVs) can span multiple genes, so we need
            # to check each gene separately, and then collapse duplicates later
            for gene in var.get_genes():
                if gene not in genes:
                    genes[gene] = []
                # add the variant to the gene entry
                genes[gene].append(var)
        
        return genes
        
    def find_variants(self, variants, gene):
        """ finds variants that fit inheritance models
        
        Args:
            variants: list of TrioGenotype objects
            gene: gene ID as string
        
        Returns:
            list of variants that pass inheritance checks
        """
        
        # get the inheritance for the gene (monoalleleic, biallelic, hemizygous
        # etc), but allow for times when we haven't specified a list of genes
        # to use
        gene_inh = None
        if self.known_genes is not None and gene in self.known_genes:
            gene_inh = self.known_genes[gene]["inh"]
        
        # If we are looking for variants in a set of known genes, and the gene
        # isn't part of that set, then we don't ant to examine the variant for
        # that gene, UNLESS the variant is a CNV, since CNVs can be included
        # purely from size thresholds, regardless of which gene they overlap.
        if self.known_genes is not None and gene not in self.known_genes:
            variants = [ x for x in variants if x.is_cnv() ]
        
        # ignore intergenic variants
        if gene is None:
            for var in variants:
                if var.get_chrom() == self.debug_chrom and var.get_position() == self.debug_pos:
                    print(var, "lacks HGNC/gene symbol")
            return []
        
        # Now that we are examining a single gene, check that the consequences
        # for the gene are in the required functional categories.
        variants = [ var for var in variants if var.child.is_lof(gene) or var.child.is_missense(gene) ]
        if variants == []:
            return []
        
        logging.debug("{} {} {} {}".format(self.family.child.get_id(), gene, variants, gene_inh))
        chrom_inheritance = variants[0].get_inheritance_type()
        
        if chrom_inheritance == "autosomal":
            finder = Autosomal(variants, self.family, self.known_genes, gene, self.cnv_regions)
        elif chrom_inheritance in ["XChrMale", "XChrFemale", "YChrMale"]:
            finder = Allosomal(variants, self.family, self.known_genes, gene, self.cnv_regions)
        
        variants = finder.get_candidate_variants()
        variants = [ (x[0], list(x[1]), list(x[2]), [gene]) for x in variants ]
        
        return variants
    
    def exclude_duplicates(self, variants):
        """ rejig variants included under multiple inheritance mechanisms
        
        Args:
            variants: list of candidate variants
        
        Returns:
            list of (variant, check_type, inheritance) tuples, with duplicates
            excluded, and originals modified to show both mechanisms
        """
        
        unique_vars = {}
        for variant in variants:
            key = variant[0].child.get_key()
            if key not in unique_vars:
                unique_vars[key] = list(variant)
            else:
                result = variant[1]
                inh = variant[2]
                hgnc = variant[3]
                
                # append the check type and inheritance type to the first
                # instance of the variant
                unique_vars[key][1] += [x for x in result if x not in unique_vars[key][1]]
                unique_vars[key][2] += [x for x in inh if x not in unique_vars[key][2]]
                
                # add the HGNC symbols that are unique to the current variant
                # to the merged variant
                hgnc = [x for x in hgnc if x not in unique_vars[key][3]]
                unique_vars[key][3] += hgnc
                
        unique_vars = [tuple(unique_vars[x]) for x in unique_vars]
        
        return unique_vars

def main():
    """ run the clinical filtering analyses
    """
    
    options = get_options()
    
    # set the level of logging to generate
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    
    log_filename = "clinical-filter.log"
    if options.ped is not None:
        log_filename = options.ped + ".log"
    
    logging.basicConfig(level=numeric_level, filename=log_filename)
    
    finder = ClinicalFilter(options)
    finder.filter_trios()

if __name__ == "__main__":
    main()
