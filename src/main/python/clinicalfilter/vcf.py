""" Loads vcf files for trios to get variants for analysis
"""

import os
import io
import platform
import sys
import gzip
import logging

from clinicalfilter import variant
from clinicalfilter.trio_genotypes import TrioGenotypes
from clinicalfilter.match_cnvs import MatchCNVs


class LoadVCFs(object):
    """ load VCF files for a trio
    """
    
    def __init__(self, pedTrio, counter, total_trios, filters):
        """ intitalise the class with the filters and trio details etc
        """
        
        self.pedTrio = pedTrio
        self.counter = counter
        self.total_trios = total_trios
        self.filters = filters
    
    def get_trio_variants(self):
        """ loads the variants for a trio
        """
        
        try:
            self.load_trio()
            self.combine_trio_variants()
            self.filter_de_novos()
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
            
            self.variants = []
        
        return self.variants
    
    def open_vcf_file(self, path):
        """Gets a VCF file handle while allowing for gzipped and text VCF formats.
        
        Args:
            path: path to VCF file.
            
        Returns:
            A file handle for VCF file.
            
        Raises:
            ValueError: An error when the VCF file path is not specified correctly.
        """
        
        if not os.path.exists(path):
            raise IOError("VCF file not found at: " + path)
        
        if path.endswith(".gz"):
            # python2 gzip opens in text, but same mode in python3 opens as bytes,
            # avoid with platform specific code
            if platform.python_version_tuple()[0] == "2":
                f = gzip.open(path, 'r')
            else:
                f = gzip.open(path, "rt")
        elif path.endswith(".vcf") or path.endswith(".txt"):
            f = io.open(path,'r', encoding="latin_1")
        else:
            print('Unable to open extension "%s". Accepted extensions are ".vcf" or ".txt"' % path.split(".")[-1])
            print('The ".txt" is tab-separated file and it can be compressed in gz format.')
            sys.exit(0)
        return f

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
    
    def split_vcf_line(self, line):
        """ splits a tab-separated vcf line
        """
        
        line = line.strip().split("\t")
        
        self.chrom = line[0]
        self.position = line[1]
        self.snp_id = line[2]
        self.ref_allele = line[3]
        self.alt_allele = line[4]
        self.quality = line[5]
        self.filter_value = line[6]
        self.info_values = line[7]
        self.format_keys = line[8]
        self.sample_values = line[9]
        
        return line
    
    def add_single_variant(self, vcf, var, gender, line):
        """ adds a single variant to a vcf dictionary indexed by position key
        """
        
        var.add_info(self.info_values)
        var.add_format(self.format_keys, self.sample_values)
        var.add_vcf_line(line)
        var.set_gender(gender)
        try:
            var.set_genotype()
            key = var.get_key()
            vcf[key] = var
        except ValueError:
            # we only get ValueError when the genotype cannot be set, which
            # occurs for x chrom male heterozygotes (an impossible genotype)
            pass
        
    
    def include_variant(self, var, child_variants, filters):
        """ check if we want to include the variant or not
        """
        
        use_variant = False
        if child_variants:
            key = var.get_key()
            if key in self.child_vcf.keys():
                use_variant = True
            elif self.alt_allele == "<DUP>" or self.alt_allele == "<DEL>":
                if self.cnv_matcher.has_match(var):
                    use_variant = True
        elif filters:
            var.add_info(self.info_values)
            if var.passes_filters(self.filters):
                use_variant = True
        else:
            use_variant = True
            
        return use_variant
        
    
    def open_individual(self, path, gender, filters=False, child_variants=False):
        """ Convert VCF to TSV format. Use for single sample VCF file.
        
        Obtains the VCF data for a single sample. This function optionally
        filters the lines of the VCF file that pass defined criteria, in order
        to reduce memory usage.
        
        Args:
            path: path to VCF file for the individual
            gender: gender of the individual
            filters: optional filters to screen variants on
            child_variants: list of tuple keys for variants identified in 
                the proband, in order to get the same variants from the 
                parents.
        
        Returns:
            A dictionary containing variant data for each variant indexed by
            a position key.
        """
        
        f = self.open_vcf_file(path)
        header = self.get_vcf_header(f)
        self.header_lines = header
        
        vcf = {}
        for line in f:
            if line.startswith("#"):
                continue
            
            line = self.split_vcf_line(line)
            
            if self.alt_allele == "<DUP>" or self.alt_allele == "<DEL>":
                var = variant.CNV(self.chrom, self.position, self.snp_id, \
                    self.ref_allele, self.alt_allele, self.quality, \
                    self.filter_value)
                var.add_info(self.info_values)
            else:
                var = variant.SNV(self.chrom, self.position, self.snp_id, \
                    self.ref_allele, self.alt_allele, self.quality, \
                    self.filter_value)
            
            # check if we want to include the variant or not
            if self.include_variant(var, child_variants, filters):
                self.add_single_variant(vcf, var, gender, line)
        
        # logging.debug(path)
        # for key in sorted(vcf):
        #     logging.debug(vcf[key])
        
        return vcf
    
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
            self.pedTrio.child.get_path())
        
        # open the childs VCF file
        self.child_vcf = self.open_individual(self.pedTrio.child.get_path(), \
            self.pedTrio.child.get_gender(), filters=True)
        self.cnv_matcher = MatchCNVs(self.child_vcf)
        
        # if the trio doesn't include parents, generate blank dictionaries
        if self.pedTrio.mother is not None:
            logging.info(" mothers path: " + self.pedTrio.mother.get_path())
            self.mother_vcf = self.open_individual(self.pedTrio.mother.get_path(), \
                self.pedTrio.mother.get_gender(), child_variants=True)
        else:
            self.mother_vcf = {}
            
        if self.pedTrio.father is not None:
            logging.info(" fathers path: " + self.pedTrio.father.get_path())
            self.father_vcf = self.open_individual(self.pedTrio.father.get_path(), \
                self.pedTrio.father.get_gender(), child_variants=True)
        else:
            self.father_vcf = {}
    
    def combine_trio_variants(self):
        """ for each variant, combine the trio's genotypes
        """
        
        mother_cnv_matcher = MatchCNVs(self.mother_vcf)
        father_cnv_matcher = MatchCNVs(self.father_vcf)
        
        self.variants = []
        for key in self.child_vcf:
            var = self.child_vcf[key]
            trio = TrioGenotypes(var)
            
            # if we only have the child, then just add the variant to the list
            if self.pedTrio.mother is None and self.pedTrio.father is None:
                self.variants.append(trio)
                continue
            
            mother_var = self.get_parental_var(var, self.mother_vcf, self.pedTrio.mother.get_gender(), mother_cnv_matcher)
            trio.add_mother_variant(mother_var)
            
            father_var = self.get_parental_var(var, self.father_vcf, self.pedTrio.father.get_gender(), father_cnv_matcher)
            trio.add_father_variant(father_var)
            
            self.variants.append(trio)
    
    def get_parental_var(self, var, parental_vcf, gender, matcher):
        """ get the corresponding parental variant to a childs variant, or 
        create a variant with reference genotype.
        
        Args:
            var: childs var, as Variant object
            parental_vcf: dict of parental variants, indexed by variant position key
            gender: gender of the parent
            matcher: cnv matcher for parent
        """
        
        if isinstance(var, variant.CNV):
            if matcher.any_overlap(var):
                overlap_key = matcher.get_overlap_key(var.get_key())
                var = parental_vcf[overlap_key]
            else:
                var = variant.CNV(var.chrom, var.position, var.id, var.ref_allele, var.alt_allele, var.quality, var.filter)
                var.set_gender(gender)
                var.set_default_genotype()
        else:
            if var.get_key() in parental_vcf:
                var = parental_vcf[var.get_key()]
            else:
                var = variant.SNV(var.chrom, var.position, var.id, var.ref_allele, var.alt_allele, var.quality, var.filter)
                var.set_gender(gender)
                var.set_default_genotype()
        
        return var
        
    def filter_de_novos(self):
        """ filter out the de novos that have been picked up at an earlier stage of the pipeline
        """
        
        # ignore situations when we haven't loaded any parents, which would all look like de novos 
        # since we insert "0" for missing parental genotypes
        if self.pedTrio.father is None and self.pedTrio.mother is None:
            return 
       
        # run through the variants in the child, and
        passed_variants = []
        for var in self.variants:
            if var.passes_de_novo_checks(self.pedTrio):
                passed_variants.append(var)
        
        self.variants = passed_variants

