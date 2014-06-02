""" Loads vcf files for trios to get variants for analysis
"""

import os
import io
import sys
import gzip
import logging
import hashlib

from clinicalfilter.variant_snv import SNV
from clinicalfilter.variant_cnv import CNV
from clinicalfilter.trio_genotypes import TrioGenotypes
from clinicalfilter.match_cnvs import MatchCNVs

IS_PYTHON2 = sys.version_info[0] == 2
IS_PYTHON3 = sys.version_info[0] == 3


class LoadVCFs(object):
    """ load VCF files for a trio
    """
    
    def __init__(self, counter, total_trios, filters, tags_dict):
        """ intitalise the class with the filters and tags details etc
        
        Args:
            counter: count of how many 
            total_trios:
            filters:
            tags_dict: 
        """
        
        self.family = None
        self.counter = counter
        self.total_trios = total_trios
        self.filters = filters
        self.tags_dict = tags_dict
    
    def get_trio_variants(self, family):
        """ loads the variants for a trio
        
        Args:
            family: Family object for a trio
        
        Returns:
            list of filtered variants for a trio, as TrioGenotypes objects
        """
        
        self.family = family
        
        try:
            (child_vars, mother_vars, father_vars) = self.load_trio()
            variants = self.combine_trio_variants(child_vars, mother_vars, father_vars)
            variants = self.filter_de_novos(variants)
        except IOError as error:
            if self.family.has_parents():
                mother_ID = self.family.mother.get_ID()
                father_ID = self.family.father.get_ID()
            else:
                mother_ID = "no mother"
                father_ID = "no father"
            
            logging.error("trio with missing file - child: " + self.family.child.get_ID() \
                + ", mother: " + mother_ID + ", father: " + father_ID + ". " + str(error))
        
        return variants
    
    def open_vcf_file(self, path):
        """ Gets a file object for an individual's VCF file.
        
        Args:
            path: path to VCF file (gzipped or text format).
            
        Returns:
            A file handle for the VCF file.
        """
        
        if not os.path.exists(path):
            raise OSError("VCF file not found at: " + path)
        
        extension = os.path.splitext(path)[1]
        
        if extension == ".gz":
            # python2 gzip opens in text, but same mode in python3 opens as
            # bytes, avoid with platform specific code
            if IS_PYTHON2:
                handle = gzip.open(path, "r")
            elif IS_PYTHON3:
                handle = gzip.open(path, "rt")
        elif extension in [".vcf", ".txt"]:
            handle = io.open(path, "r", encoding="latin_1")
        else:
            raise OSError("unsupported filetype: " + path)
        
        return handle

    def get_vcf_header(self, f):
        """Get the header lines from a VCF file.
        
        Args:
            f: file object to read from
        
        Returns:
            a list of lines that start with "#", which are the header lines.
        """
        header = []
        
        for line in f:
            if not line.startswith("#"):
                break
            header.append(line)
        
        return header
    
    def add_single_variant(self, variants, var, gender, line):
        """ adds a single variant to a vcf dictionary indexed by position key
        
        Args:
            variants: list of variants for an individual
            var: single Variant object
            gender: gender of the individual (eg "M"/"F" or "1"/"2")
            line: list of elements of the VCF line for the variant
        """
        
        var.add_format(line[8], line[9])
        var.add_vcf_line(line)
        var.set_gender(gender)
        
        try:
            var.set_genotype()
            variants.append(var)
        except ValueError:
            # we only get ValueError when the genotype cannot be set, which
            # occurs for x chrom male heterozygotes (an impossible genotype)
            pass
    
    def find_vcf_definitions(self, path, header):
        """ get provenance information for a vcf path
        
        Args:
            path: path to VCF file
            header: list of header lines
        
        Returns:
            returns a tuple of sha1 VCF file hash, name of VCF file (without 
            directory), and date the VCF file was generated
        """ 
        
        try:
            vcf_checksum = hashlib.sha1(open(path, "rb").read()).hexdigest()
        except (OSError, IOError) as e:
            vcf_checksum = "NA"
        
        vcf_basename = os.path.basename(path)
        
        vcf_date = None
        for line in header:
            if line.startswith("##fileDate"):
                vcf_date = line.strip().split("=")[1]
                break
        
        if vcf_date is None:
            vcf_date = vcf_basename.split(".")[-3]
        
        return (vcf_checksum, vcf_basename, vcf_date)
    
    def construct_variant(self, line, gender):
        """ constructs a Variant object for a VCF line, specific to the variant type
        
        Args:
            line: list of elements of a single sample VCF line:
                [chrom, position, snp_id, ref_allele, alt_allele, quality, 
                filter_value, info, format_keys, format_values]
            gender: gender of the individual to whom the variant line belongs 
                (eg "1" or "M" for male, "2", or "F" for female).
        
        Returns:
            returns a Variant object
        """
        
        # CNVs are found by their alt_allele values, as either <DUP>, or <DEL>
        if line[4] == "<DUP>" or line[4] == "<DEL>":
            var = CNV(line[0], line[1], line[2], line[3], line[4], line[6])
            var.add_info(line[7], self.tags_dict)
            # CNVs require the format values for filtering
            var.set_gender(gender)
            var.add_format(line[8], line[9])
            if "HGNC" in self.filters:
                var.fix_gene_IDs(self.filters["HGNC"][1])
        else:
            var = SNV(line[0], line[1], line[2], line[3], line[4], line[6])
            var.add_info(line[7], self.tags_dict)
        
        return var
    
    def include_variant(self, line, child_variants, gender):
        """ check if we want to include the variant or not
        
        Args:
            line: list of elements from the VCF line for the variant.
            child_variants: True/False for whether variants have been filtered
                for the proband (if so, we can simply check the parent's 
                variants for matches in the child's variants).
            gender: the gender of the proband (used in CNV filtering).
        
        Returns:
            True/False for whether to include the variant.
        """
        
        use_variant = False
        if child_variants:
            key = (line[0], line[1])
            if key in self.child_vcf:
                use_variant = True
            elif line[4] == "<DUP>" or line[4] == "<DEL>":
                var = self.construct_variant(line, gender)
                if self.cnv_matcher.has_match(var):
                    use_variant = True
        else:
            var = self.construct_variant(line, gender)
            if var.passes_filters(self.filters):
                use_variant = True
            
        return use_variant
        
    def open_individual(self, individual, child_variants=False):
        """ Convert VCF to TSV format. Use for single sample VCF file.
        
        Obtains the VCF data for a single sample. This function optionally
        filters the lines of the VCF file that pass defined criteria, in order
        to reduce memory usage.
        
        Args:
            individual: Person object for individual
            child_variants: True/False for whether variants have been filtered
                for the proband (if so, we can simply check the parent's 
                variants for matches in the child's variants).
        
        Returns:
            A dictionary containing variant data for each variant indexed by
            a position key.
        """
        
        path = individual.get_path()
        gender = individual.get_gender()
        
        vcf = self.open_vcf_file(path)
        self.header_lines = self.get_vcf_header(vcf)
        file_definitions = self.find_vcf_definitions(path, self.header_lines)
        
        variants = []
        for line in vcf:
            line = line.strip().split("\t")
            
            # check if we want to include the variant or not
            if self.include_variant(line, child_variants, gender):
                var = self.construct_variant(line, gender)
                self.add_single_variant(variants, var, gender, line)
        
        return variants, file_definitions
    
    def load_trio(self):
        """ opens and parses the VCF files for members of the family trio.
        
        We need to load the VCF data for each of the members of the trio. As a
        bare minimum we need VCF data for the child in the family. Occasionally
        we lack parents for the child, so we create blank entries when that 
        happens.
        """
        
        # load the VCF file for each member of the trio
        logging.info("opening trio " + str(self.counter + 1) + " of " + \
            str(self.total_trios) + ". child path: " + \
            self.family.child.get_path())
        
        # open the childs VCF file
        child_vars, self.child_defs = self.open_individual(self.family.child)
        
        # get the childs variant keys, so we can check if they are in the 
        # parents VCF
        self.child_vcf = set()
        for var in child_vars:
            self.child_vcf.add(var.get_key())
        
        self.child_header = self.header_lines
        self.cnv_matcher = MatchCNVs(child_vars)
        
        if self.family.has_parents():
            logging.info(" mothers path: " + self.family.mother.get_path())
            mother_vars, self.mother_defs = self.open_individual(self.family.mother, child_variants=True)
            
            logging.info(" fathers path: " + self.family.father.get_path())
            father_vars, self.father_defs = self.open_individual(self.family.father, child_variants=True)
        else:
            # if the trio doesn't include parents, generate blank dictionaries
            mother_vars, self.mother_defs = [], ("NA", "NA", "NA")
            father_vars, self.father_defs = [], ("NA", "NA", "NA")
        
        return (child_vars, mother_vars, father_vars)
    
    def combine_trio_variants(self, child_vars, mother_vars, father_vars):
        """ for each variant, combine the trio's genotypes
        """
        
        mother_cnv_matcher = MatchCNVs(mother_vars)
        father_cnv_matcher = MatchCNVs(father_vars)
        
        variants = []
        for var in child_vars:
            trio = TrioGenotypes(var)
            
            # if we only have the child, then just add the variant to the list
            if self.family.has_parents() == False:
                variants.append(trio)
                continue
            
            mother_var = self.get_parental_var(var, mother_vars, self.family.mother.get_gender(), mother_cnv_matcher)
            trio.add_mother_variant(mother_var)
            
            father_var = self.get_parental_var(var, father_vars, self.family.father.get_gender(), father_cnv_matcher)
            trio.add_father_variant(father_var)
            
            variants.append(trio)
        
        return variants
    
    def get_parental_var(self, var, parental_vcf, gender, matcher):
        """ get the corresponding parental variant to a childs variant, or 
        create a default variant with reference genotype.
        
        Args:
            var: childs var, as Variant object
            parental_vcf: dict of parental variants, indexed by position key
            gender: gender of the parent
            matcher: cnv matcher for parent
        
        Returns:
            returns a Variant object, matched to the proband's variant
        """
        
        var_key = var.get_key()
        
        # if the variant is a CNV, the corresponding variant might not match 
        # the start site, so we look a variant that overlaps
        if isinstance(var, CNV) and matcher.has_match(var):
            var_key = matcher.get_overlap_key(var.get_key())
            
        for parental_var in parental_vcf:
            if var_key == parental_var.get_key():
                return parental_var
        
        # if the childs variant does not exist in the parents VCF, then we 
        # create a default variant for the parent
        if isinstance(var, CNV):
            var = CNV(var.chrom, var.position, var.id, var.ref_allele, var.alt_allele, var.filter)
        else:
            var = SNV(var.chrom, var.position, var.id, var.ref_allele, var.alt_allele, var.filter)
        
        var.set_gender(gender)
        var.set_default_genotype()
        
        return var
    
    def get_vcf_provenance(self):
        """ returns provenance of VCFs for individuals in a trio
        """
        
        return self.child_defs, self.mother_defs, self.father_defs
    
    def filter_de_novos(self, variants):
        """ filter the de novos variants in the VCF files
        
        Args:
            variants: list of TrioGenotypes objects.
        
        Returns:
            a list of TrioGenotypes without the de novo variants that failed the 
            de novo filter.
        """
        
        # ignore situations when we haven't loaded any parents, which would all
        # look like de novos, since we insert "0" for missing parental genotypes
        if self.family.has_parents() == False:
            return 
       
        # run through the variants in the child, and remove de novos that fail
        # denovogear filtering criteria
        passed_variants = []
        for var in variants:
            if var.passes_de_novo_checks(self.family):
                passed_variants.append(var)
        
        return passed_variants

