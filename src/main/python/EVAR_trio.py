''' Clinical filtering for trios

Find variants in affected children that might contribute to their disorder. We
load VCF files (either named on the command line, or listed in a PED file) for
members of a family, filter for rare, functionally disruptive variants, and
assess whether each variant might affect the child's disorder. We take into 
account the parents genotypes (if available) and whether the parents are also
affected with a (the?) disorder. For variants in known disease causative genes
we check whether the inheritance patterns matches one expected for the  
inheritance models of the gene.

Usage:

python EVAR_trio.py \
    --ped temp_name.ped \
    --filter filters.txt \
    --tags tags.txt \
    --known-genes known_genes.txt \
    --alternate-ids alternate_ids.txt \
    --output output_name.txt


sa9@sanger.ac.uk
'''

import platform
if platform.python_implementation() == "PyPy":
    import numpypy

import numpy as np
import sys
import optparse
import itertools
import os
import logging

from clinicalfilter import user
from clinicalfilter import vcf
from clinicalfilter import parser
from clinicalfilter import inheritance
from clinicalfilter import ped
from clinicalfilter import reporting

pseudoautosomal_regions = [(1,2699520), (154930290,155260560), (88456802,92375509)]

class Trio(parser.Parser, reporting.report):
    def __init__(self, defs):
        
        # set some definitions
        self.setDefinitions(defs)
        self.first_run = True
        
        self.counter = 0
        # load the trio paths into the current path setup, currently only test a single trio
        for family_ID in sorted(self.pedTrios):
            self.pedTrio = self.pedTrios[family_ID]
            
            # some families have more than one child in the family, so run through each child.
            while self.pedTrio.check_all_children_analysed() == False:
                self.pedTrio.set_child()
                if self.pedTrio.child.get_boolean_affected_status():
                    try:
                        self.load_vcfs()
                        self.filter_de_novos()
                        self.analyse_trio()
                    except IOError as error:
                        if self.pedTrio.mother is None:
                            mother_ID = "no mother"
                        else:
                            mother_ID = self.pedTrio.mother.get_ID()
                        if self.pedTrio.father is None:
                            father_ID = "no father"
                        else:
                            father_ID = self.pedTrio.father.get_ID()
                        logging.error("trio with missing file - child: " + self.pedTrio.child.get_ID() \
                            + ", mother: " + mother_ID + ", father: " + father_ID + ". " + str(error))
                
                self.pedTrio.set_child_examined()
            self.counter += 1
        
        # make sure we close the output file off
        if self.output_path is not None:
            self.output.close()
        
        sys.exit(0)
    
    def analyse_trio(self):
        """identify candidate variants in exome data for a single trio.
        """
        
        # report the variants that were found
        if self.first_run:
            self.printReportInputInfo('\tTrio family analysis [EVA Report]')
            
            # open the output file while analysing the first trio
            if self.output_path is not None:
                self.output = open(self.output_path, "w")
        
        # create numpy matrices of variants and genotypes in each gene, then find variants that fit
        # different inheritance models
        self.create_gene_matrix()
        self.found_variants = {}
        for gene in self.genes_dict:
            chrom = self.get_chr(self.genes_dict[gene]["positions"])
            positions = self.get_sorted_positions(self.genes_dict[gene]["positions"].keys())
            chr_inheritance_group = self.get_chr_group(chrom, positions[0], positions[-1])
            
            # get the inheritance for the gene (monoalleleic, biallelic, hemizygous etc), but allow
            # for times when we haven't specified a list of genes to use
            if self.known_genes is not None:
                gene_inheritance = self.known_genes[gene]
            else:
                gene_inheritance = None
            
            # ignore intergenic variants
            if gene == None:
                continue
            
            logging.debug(self.pedTrio.child.get_ID() + " " + gene + " " + str(self.genes_dict[gene]["positions"].keys()) + " " + str(gene_inheritance))
            finder = inheritance.inheritance(self.genes_dict[gene]["positions"], self.pedTrio, chr_inheritance_group, gene_inheritance)
            candidates = finder.get_candidiate_variants()
            for candidate in candidates:
                self.include_variant(gene, candidate)
        
        # export the results to either tab-separated table or VCF format
        if self.output_path is not None:
            self.save_results()
        if self.export_vcf:
            self.save_vcf()
        self.first_run = False
    
    def load_vcfs(self):
        """ opens and parses the VCF files for members of the family trio.
        
        We need to load the VCF data for each of the members of the trio. As a bare minimum we need 
        VCF data for the child in the family. Occasionally we lack parents for the child, so we 
        create blank entries when that happens.
        """
        
        # load the VCF file for each member of the trio
        logging.info("opening trio " + str(self.counter + 1) + " of " + str(len(self.pedTrios)) \
                     + ". child path: " + self.pedTrio.child.get_path())
        
        # open the childs VCF file
        self.child_vcf = vcf.vcf2tsv(self.pedTrio.child.get_path(), self.isPassUserFilters)
        
        # if the child doesn't have a parent listed, generate blank dictionary for them
        if self.pedTrio.mother is not None:
            logging.info(" mothers path: " + self.pedTrio.mother.get_path())
            self.mother_vcf = vcf.vcf2tsv(self.pedTrio.mother.get_path(), child_variants=self.child_vcf["data"].keys())
        else:
            self.mother_vcf = {"data": {}}
            
        if self.pedTrio.father is not None:
            logging.info(" fathers path: " + self.pedTrio.father.get_path())
            self.father_vcf = vcf.vcf2tsv(self.pedTrio.father.get_path(), child_variants=self.child_vcf["data"].keys())
        else:
            self.father_vcf = {"data": {}}
    
    def setDefinitions(self, defs):
        """sets the definitions for the filters and weights, as previously created by the 
        loadDefinitions class.
        """
        
        self.filters_path = defs.filters_path
        self.output_path = defs.output_path
        self.export_vcf = defs.export_vcf
        
        self.filters = defs.filters
        self.tags_dict = defs.tags_dict
        
        self.known_genes = defs.known_genes
        self.ID_mapper = defs.ID_mapper
        
        self.pedTrios = defs.pedTrios
    
    def create_gene_matrix(self):
        """creates numpy arrays of trio variant genotypes and weights for each gene.
        
        For each gene in self.genes_dict, extract the genotypes for the variants, then convert 
        the genotypes for a variant into a numpy array, along with the weights. This seems like a 
        overly complex way of going about it.
        """
        
        # organise the variants into entries for each gene
        self.genes_dict = {}
        self.make_genes_dict(self.child_vcf, "child")
        self.make_genes_dict(self.mother_vcf, "mother")
        self.make_genes_dict(self.father_vcf, "father")
        
        for gene in self.genes_dict:
            sorted_pos = self.get_sorted_positions(self.genes_dict[gene]["positions"])

            family_members = ['child', 'mother', 'father']
            for pos in sorted_pos:
                for member in family_members:
                    # extract the genotype and weights for the variant, unless the individual does 
                    # not have them
                    if member in self.genes_dict[gene]["positions"][pos]:
                        record = self.genes_dict[gene]["positions"][pos][member]
                        # get the genotype, but make sure we account for different ways of naming
                        # the genotype field
                        genotype = None
                        for genotype_tag in self.tags_dict["genotype"]:
                            if genotype_tag in record:
                                try:
                                    genotype = vcf.translateGT(record[genotype_tag])
                                except TypeError:
                                    genotype = None
                                break
                    else:
                        # the member does not have variants in this pos i.e. hom ref
                        genotype = 0
                        self.genes_dict[gene]["positions"][pos][member] = {}
                    
                    self.genes_dict[gene]["positions"][pos][member]["genotype"] = genotype
    
    def inPARRegion(self, pos):
        """checks whether a nucleotide position lies within defined pseudoautosomal regions.
        
        Args:
            pos: nucleotide position in bp_position
        
        Returns:
            True false for whether the position falls within a pseudoautosomal region
        """
        pos = int(pos)
        for start, end in pseudoautosomal_regions:
            if pos >= start and pos <= end:
                return True
        return False
    
    def get_chr(self, record):
        """returns the chromosome ID for a gene.
        
        Args:
            record: dictionary entry for a gene, containing all the variants
            
        Returns:
            the chromosome that the gene is on.
        """
        for pos in record:
            for member in ['child', 'mother', 'father']:
                if member in record[pos]:
                    chrom = record[pos][member]['CHROM']
                    break
        return chrom

    def get_chr_group(self, chrom, min_pos, max_pos):
        """ returns the chromosome type (eg autosomal, or X chromosome type).
        
        provides the chromosome type for a chromosome (eg Autosomal, or X-chrom male etc). This only
        does simple string matching. The chromosome string is either the chromosome number, or in 
        the case of the sex-chromosomes, the chromosome character. This doesn't allow for 
        chromosomes to be specified as 'chr1', and sex chromosomes have to be specified as 'X' or 
        'Y', not '23' or '24'.
        
        Args:
            CHR: chromosome string, coded just as the ID (eg '1', '2', '3' ... 'X', 'Y')
        
        Returns:
            chromosome type (eg 'autosomal', 'XChrMale', 'XChrFemale')
        """
        
        if chrom not in ['chrX', 'ChrX', 'X']:
            return 'autosomal'
        else:
            # check if the gene lies within a pseudoautosomal region
            for start, end in pseudoautosomal_regions:
                if start < int(min_pos) < end or start < int(max_pos) < end:
                    return "autosomal"
            
            if self.pedTrio.child.is_male():
                return 'XChrMale'
            if self.pedTrio.child.is_female():
                return 'XChrFemale'  
    
    def get_sorted_positions(self, positions):
        """Returns a list of nucleotide positions, sorted according to the nucleotide position.
        
        Args:
            positions: list of nucleotide position strings
        
        Returns:
            list of sorted nucleotide positions, as strings.
        """
        
        sorted_pos = sorted([int(pos) for pos in positions])
        sorted_pos = list(map(str, sorted_pos))
        
        return sorted_pos     
    
    def include_variant(self, gene, candidate):
        """Adds a candidate variant to dictionary of variants for a specific inheritance model.
        
        Args:
            gene: HUGO gene ID
            candidate: a list of [variant, nucleotide position, check_type, inheritance]
        """
        
        dictionary = self.found_variants
        
        variant = candidate[0]
        position = candidate[1]
        check_type = candidate[2]
        inheritance_type = candidate[3]
        record = self.extractUserFields(variant["child"])
        record["MAX_MAF"] = self.find_max_allele_frequency(variant["child"])
        
        variant["child"]["genotype"] = vcf.translateGT(variant["child"]["genotype"])
        
        trio_genotype = "%d/%d/%d" % (variant["child"]["genotype"], variant["mother"]["genotype"], variant["father"]["genotype"])
        record['trio_genotype'] = trio_genotype
        
        # make sure the gene is in the dictionary
        if gene not in dictionary:
            dictionary[gene] = {}
        
        # make sure the VCF record for the child is recorded under the variant position entry for 
        # the gene
        dictionary[gene][position] = record
        
        # add in all sorts of extra information here, probably better to do this later though
        dictionary[gene][position]['inh'] = inheritance_type
        
        if self.pedTrio.father is not None:
            dictionary[gene][position]['dad_aff'] = self.pedTrio.father.get_affected_status()
        else:
            dictionary[gene][position]['dad_aff'] = "NA"
        if self.pedTrio.mother is not None:
            dictionary[gene][position]['mom_aff'] = self.pedTrio.mother.get_affected_status()
        else:
            dictionary[gene][position]['mom_aff'] = "NA"
        
        dictionary[gene][position]['person_ID'] = self.pedTrio.child.get_ID()
        dictionary[gene][position]['result'] = check_type
        
        # include an alternate ID for the affected child, if it exists
        if self.ID_mapper is not None:
            dictionary[gene][position]['alternate_ID'] = self.ID_mapper[self.pedTrio.child.get_ID()]
        else:
            dictionary[gene][position]['alternate_ID'] = 'no_alternate_ID'
    
    def get_genotype(self, record):
        """extracts the genotype for an individual for a VCF record
        
        Args:
            record: a VCF record as a dictionary for a single variant for a single person
        
        Returns:
            A genotype converted to 0, 1, or 2 (referring to number of copies of the alternate 
            allele). If the record does not have a genotype, then return None.
        """
        
        for genotype_tag in self.tags_dict["genotype"]:
            if genotype_tag in record:
                try:
                    return vcf.translateGT(record[genotype_tag])
                except TypeError:
                    return None
        else:
            return None
    
    def filter_de_novos(self):
        """ filter out the de novos that have been picked up at an earlier stage of the pipeline
        """
        
        # ignore situations when we haven't loaded any parents, which would all look like de novos 
        # since we insert "0" for missing parental genotypes
        if self.pedTrio.father is None and self.pedTrio.mother is None:
            return
       
        # run through the variants in the child, and
        positions = list(self.child_vcf["data"].keys())
        for position in positions:
            if self.filter_de_novo_for_single_variant(position, self.child_vcf, self.mother_vcf, self.father_vcf) == False:
                del self.child_vcf["data"][position]
                if position in self.mother_vcf["data"]:
                    del self.mother_vcf["data"][position]
                if position in self.father_vcf["data"]:
                    del self.father_vcf["data"][position]
    
    def filter_de_novo_for_single_variant(self, position, child_vcf, mother_vcf, father_vcf):
        """ check if a variant is de novo in the child, and whether screened by denovogear
        
        Some variants are de novo in the child, and if they are, then we should subject them to 
        additional filtering to see if they have been passed by de novo gear, and an additional hard
        coded filter called 'TEAM29_FILTER' that describes whether the variant passed screening, or 
        if not, which filter it failed.
        
        Args:
            position: chromosome and nucleotide position tuple for the variant
            child_vcf: full vcf data for the child
            mother_vcf: full vcf data for the mother
            father_vcf: full vcf data for the father
        
        Returns:
            boolean value for whether the variant should be included
        """
        
        # currently hard code the filtering fields. The de novo field indicates whether the variant
        # is de novo, I don't know how this is assigned. The project filter field indicates an 
        # internal filter, curently whether the variant passed MAF, alternate frequency, and 
        # segmental duplication criteria.
        de_novo_snp_field = "DENOVO-SNP"
        de_novo_indel_field = "DENOVO-INDEL"
        project_filter_field = "TEAM29_FILTER"
        
        # set the standard de novo genotype combination
        de_novo_genotype = (1,0,0)
        
        # account for X chrom de novos in males
        chrom = position[0]
        bp_position = position[1]
        chrom_inheritance = self.get_chr_group(chrom, bp_position, bp_position)
        if chrom_inheritance == "XChrMale" and self.pedTrio.child.is_male():
            de_novo_genotype = (2,0,0)
        
        # get the genotypes for the trio
        trio = [child_vcf, mother_vcf, father_vcf]
        trio_genotype = []
        for member in trio:
            if position in member["data"]:
                genotype = self.get_genotype(member["data"][position])
                if genotype is None:
                    genotype = 0
                trio_genotype.append(genotype)
            else:
                trio_genotype.append(0)
        
        trio_genotype = tuple(trio_genotype)
        
        # if the variant is not de novo, don't worry about the additional filtering
        if trio_genotype != de_novo_genotype:
            return True
        
        # check the childs VCF record to see whether the variant has been screened out
        child_vcf = child_vcf["data"][position]
        if child_vcf[de_novo_snp_field] is None and child_vcf[de_novo_indel_field] is None:
            return False
        
        if child_vcf[project_filter_field] != "PASS":
            return False
        else:
            return True


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
        self.filters = user.parseFilters(self.filters_path)
        self.tags_dict = user.parseTags(self.tags_path)
        
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
            self.known_genes = user.parseGenes(self.options.genes_path)
            # include all the possible ways IDs that a gene field can be named in a VCF file
            for tag in self.tags_dict['gene']:
                self.filters[tag] = ["list", self.known_genes]
        else:
            self.known_genes = None
        
        # if we have named an ID mapping file, the load a dictionary of IDs and alternate IDs, so we
        # can convert between different ID schemes.
        if self.options.alternate_ids_path is not None:
            self.ID_mapper = user.create_person_ID_mapper(self.options.alternate_ids_path)
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
            pedTrio = ped.pedTrio('blank_family_ID')
            pedTrio.set_child('child', self.options.child_path, "2", self.options.child_gender)
            if self.options.mother_path is not None:
                pedTrio.set_mother('mother', self.options.mother_path, self.options.mother_affected, '2')
            if self.options.father_path is not None:
                pedTrio.set_father('father', self.options.father_path, self.options.father_affected, '1')
            
            self.pedTrios = [pedTrio]
        else:
            self.pedTrios = ped.loadPedTrios(self.options.ped_path)


def get_options():
    """gets the options from the command line
    """
    
    parser = optparse.OptionParser()
    parser.add_option('-p', '--ped', dest='ped_path', help='path to ped file containing cohort details for multiple trios')
    parser.add_option('-c', '--child', dest='child_path', help='path to child\'s VCF file')
    parser.add_option('-m', '--mother', dest='mother_path', help='path to mother\'s VCF file')
    parser.add_option('-f', '--father', dest='father_path', help='path to father\'s VCF file')
    parser.add_option('-G', '--gender', dest='child_gender', help='The child gender (male or female)')
    parser.add_option('--mom_aff_status', dest='mother_affected', help='affected status of the mother (1 = unaffacted, or 2 = affected)')
    parser.add_option('--dad_aff_status', dest='father_affected', help='affected status of the father (1 = unaffacted, or 2 = affected)')
    
    parser.add_option('-l', '--filter', dest='filters_path', help='path to filter file (eg filters.txt)')
    parser.add_option('-t', '--tags', dest='tags_path', help='path to tags.txt (eg tags.txt)')
    parser.add_option('--known-genes', dest='genes_path', help='path to list of known disease causative genes, eg DDG2P-reportable.txt')
    parser.add_option('--alternate-ids', dest='alternate_ids_path', help='path to list of alternate IDs, eg personid_decipher_id_sangerid.txt')
    parser.add_option('-o', '--output', dest='output_path', default='clinical_reporting.txt', help='filename to output variant data to')
    parser.add_option('--export-vcf', dest='export_vcf', default=False, help='whether to export identified variants to a VCF file')
    parser.add_option('--log', dest='loglevel', default="debug", help='level of logging to use, choose from: debug, info, warning, error or critical')
    
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
        raise ValueError('Invalid log level: %s' % loglevel)
    if opts.ped_path is not None:
        log_filename = opts.ped_path + '.log'
    else:
        log_filename = "clinical-filter.log"
    logging.basicConfig(level=numeric_level, filename=log_filename)
    
    defs = loadDefinitions(opts)
    myTrio = Trio(defs)

if __name__ == "__main__":
    main()

