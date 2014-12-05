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
    
    def __init__(self, total_trios, filters, tags_dict, debug_chrom, debug_pos):
        """ intitalise the class with the filters and tags details etc
        
        Args:
            counter: count of how many trios have been analysed
            total_trios: count of how many trios are to be analysed
            filters: dictionary of filtering criteria for variants
            tags_dict: dictionary of alternate tags for INFO fields
        """
        
        self.family = None
        self.counter = 0
        self.total_trios = total_trios
        self.filters = filters
        self.tags_dict = tags_dict
        
        SNV.debug_chrom = debug_chrom
        SNV.debug_pos = debug_pos
        CNV.debug_chrom = debug_chrom
        CNV.debug_pos = debug_pos
        
        if debug_chrom is not None:
            SNV.passes_filters = SNV.passes_filters_with_debug
    
    def get_trio_variants(self, family, pp_filter):
        """ loads the variants for a trio
        
        Args:
            family: Family object for a trio
            pp_filter float between 0 and 1, being the threshold for the PP_DNM filter
        
        Returns:
            list of filtered variants for a trio, as TrioGenotypes objects
        """
        
        self.family = family
        self.counter += 1
        
        try:
            (child_vars, mother_vars, father_vars) = self.load_trio()
            variants = self.combine_trio_variants(child_vars, mother_vars, father_vars)
            variants = self.filter_de_novos(variants, pp_filter)
        except IOError as error:
            if self.family.has_parents():
                mother_id = self.family.mother.get_id()
                father_id = self.family.father.get_id()
            else:
                mother_id = "no mother"
                father_id = "no father"
            
            logging.error("trio with missing file - child: " + self.family.child.get_id() \
                + ", mother: " + mother_id + ", father: " + father_id + ". " + str(error))
            
            raise(error)
        
        print(variants)
        
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

    def get_vcf_header(self, path):
        """ Get the header lines from a VCF file.
        
        Args:
            path: path to VCF file
        
        Returns:
            a list of lines that start with "#", which are the header lines.
        """
        
        f = self.open_vcf_file(path)
        
        header = []
        for line in f:
            if not line.startswith("#"):
                break
            header.append(line)
        
        f.close()
        
        return header
    
    def exclude_header(self, f):
        """ removes the header from a VCF file object
        
        We remove the header from the VCF file, since the header is ~200 lines
        long, and an exome VCF file is 100,000 lines long, so it's better to 
        remove the header once, rather than continually check if lines are part 
        of the header as we traverse the VCF. We simply run through the VCF 
        until we find a non-header line, then seek back to the start of that 
        line.
        
        Args:
            f: handler for a VCF file
        """
        
        file_pos = f.tell()
        in_header = True
        while in_header:
            line = f.readline()
            if not line.startswith("#"):
                in_header = False
                break
            file_pos = f.tell()   
        
        # jump back to the start of the last line, which should be the first
        # line following the header
        f.seek(file_pos)
    
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
            if key in self.child_keys:
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
            A list of variants for the individual.
        """
        
        path = individual.get_path()
        gender = individual.get_gender()
        
        # open the vcf, and adjust the position in the file to immediately after 
        # the header, so we can run through the variants
        vcf = self.open_vcf_file(path)
        self.exclude_header(vcf)
        
        variants = []
        for line in vcf:
            line = line.strip().split("\t")
            
            # check if we want to include the variant or not
            if self.include_variant(line, child_variants, gender):
                var = self.construct_variant(line, gender)
                self.add_single_variant(variants, var, gender, line)
        
        return variants
    
    def load_trio(self):
        """ opens and parses the VCF files for members of the family trio.
        
        We need to load the VCF data for each of the members of the trio. As a
        bare minimum we need VCF data for the child in the family. Occasionally
        we lack parents for the child, so we create blank entries when that 
        happens.
        """
        
        # load the VCF file for each member of the trio
        logging.info("opening trio " + str(self.counter) + " of " + \
            str(self.total_trios) + ". child path: " + \
            self.family.child.get_path())
        
        # open the childs VCF file, and get the variant keys, to check if they
        # are in the parents VCF
        child_vars = self.open_individual(self.family.child)
        self.child_keys = set([var.get_key() for var in child_vars])
        
        self.child_header = self.get_vcf_header(self.family.child.get_path())
        self.cnv_matcher = MatchCNVs(child_vars)
        
        mother_vars = []
        father_vars = []
        if self.family.has_parents():
            logging.info(" mothers path: " + self.family.mother.get_path())
            mother_vars = self.open_individual(self.family.mother, child_variants=True)
            
            logging.info(" fathers path: " + self.family.father.get_path())
            father_vars = self.open_individual(self.family.father, child_variants=True)
        
        return (child_vars, mother_vars, father_vars)
    
    def combine_trio_variants(self, child_vars, mother_vars, father_vars):
        """ for each variant, combine the trio's genotypes into TrioGenotypes
        
        Args:
            child_vars: list of Variant objects for the child
            mother_vars: list of Variant objects for the mother
            father_vars: list of Variant objects for the father
        
        Returns:
            list of TrioGenotypes objects for the family
        """
        
        mother_cnv_matcher = MatchCNVs(mother_vars)
        father_cnv_matcher = MatchCNVs(father_vars)
        
        variants = []
        for var in child_vars:
            trio = TrioGenotypes(var, SNV.debug_chrom, SNV.debug_pos)
            
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
    
    def get_parental_var(self, var, parental_vars, gender, matcher):
        """ get the corresponding parental variant to a childs variant, or 
        create a default variant with reference genotype.
        
        Args:
            var: childs var, as Variant object
            parental_vars: list of parental variants
            gender: gender of the parent
            matcher: cnv matcher for parent
        
        Returns:
            returns a Variant object, matched to the proband's variant
        """
        
        key = var.get_key()
        
        # if the variant is a CNV, the corresponding variant might not match 
        # the start site, so we look a variant that overlaps
        if isinstance(var, CNV) and matcher.has_match(var):
            key = matcher.get_overlap_key(key)
            
        for parental in parental_vars:
            if key == parental.get_key():
                return parental
        
        # if the childs variant does not exist in the parents VCF, then we 
        # create a default variant for the parent
        if isinstance(var, CNV):
            parental = CNV(var.chrom, var.position, var.id, var.ref_allele, var.alt_allele, var.filter)
        else:
            parental = SNV(var.chrom, var.position, var.id, var.ref_allele, var.alt_allele, var.filter)
        
        parental.set_gender(gender)
        parental.set_default_genotype()
        
        return parental
    
    def get_vcf_provenance(self, path):
        """ get provenance information for a vcf path
        
        Args:
            path: path to VCF file
        
        Returns:
            returns a tuple of sha1 VCF file hash, name of VCF file (without
            directory), and date the VCF file was generated
        """
        
        handle = open(path, "rb")
        vcf_checksum = hashlib.sha1(handle.read()).hexdigest()
        handle.close()
        
        vcf_basename = os.path.basename(path)
        
        header = self.get_vcf_header(path)
        
        vcf_date = None
        for line in header:
            if line.startswith("##fileDate"):
                vcf_date = line.strip().split("=")[1]
                break
        
        # some VCF files lack the fileDate in the header, get it from the path
        if vcf_date is None:
            vcf_date = os.path.splitext(vcf_basename)[0]
            vcf_date = vcf_date.split(".")[2]
        
        return (vcf_checksum, vcf_basename, vcf_date)
    
    def get_trio_provenance(self):
        """ returns provenance of VCFs for individuals in a trio
        """
        
        child_defs = self.get_vcf_provenance(self.family.child.get_path())
        
        mother_defs = ("NA", "NA", "NA")
        father_defs = ("NA", "NA", "NA")
        if self.family.has_parents():
            mother_defs = self.get_vcf_provenance(self.family.mother.get_path())
            father_defs = self.get_vcf_provenance(self.family.father.get_path())
        
        return child_defs, mother_defs, father_defs
    
    def filter_de_novos(self, variants, pp_filter):
        """ filter the de novos variants in the VCF files
        
        Args:
            variants: list of TrioGenotypes objects.
            pp_filter float between 0 and 1, being the threshold for the PP_DNM filter
        
        Returns:
            a list of TrioGenotypes without the de novo variants that failed the 
            de novo filter.
        """
        
        # ignore situations when we haven't loaded any parents, which would all
        # look like de novos, since we insert "0" for missing parental genotypes
        if self.family.has_parents() == False:
            return variants
       
        # run through the variants in the child, and remove de novos that fail
        # denovogear filtering criteria
        passed_variants = []
        for var in variants:
            if var.passes_de_novo_checks(pp_filter):
                passed_variants.append(var)
        
        return passed_variants

