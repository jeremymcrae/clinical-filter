'''
Copyright (c) 2016 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import logging

from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.trio_genotypes import TrioGenotypes
from clinicalfilter.match_cnvs import MatchCNVs
from clinicalfilter.utils import open_vcf, get_vcf_header, exclude_header, \
    construct_variant
from clinicalfilter.multinucleotide_variants import get_mnv_candidates

class LoadVCFs(object):
    """ load VCF files for a trio
    """
    
    def __init__(self, total_trios, known_genes, last_base, debug_chrom, debug_pos):
        """ intitalise the class with the filters and tags details etc
        
        Args:
            total_trios: count of how many trios are to be analysed
            known_genes: dictionary of genes known to be involved with genetic
                disorders.
            last_base: set of sites in genome at conserved last base of exons,
                where we upgrade the severity of variants to loss-of-function.
            debug_chrom: chromosome string, to give more information about why
                a variant fails to pass the filters.
            debug_pos: chromosome position, to give more information about why
                a variant fails to pass the filters.
        """
        
        self.family = None
        self.counter = 0
        self.total_trios = total_trios
        self.known_genes = known_genes
        
        # define several parameters of the variant classes, before we have
        # initialised any class objects
        SNV.set_known_genes(known_genes)
        SNV.set_debug(debug_chrom, debug_pos)
        SNV.set_last_base_sites(last_base)
        
        CNV.set_known_genes(known_genes)
        CNV.set_debug(debug_chrom, debug_pos)
    
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
        logging.info("opening trio {} of {}".format(self.counter, self.total_trios))
        
        try:
            variants = self.load_trio(family)
            variants = self.filter_de_novos(variants, pp_filter)
        except OSError as error:
            mother_id, father_id = "no mother", "no father"
            if family.has_parents():
                mother_id = family.mother.get_id()
                father_id = family.father.get_id()
            
            logging.error("trio with missing file - child: " + family.child.get_id() \
                + ", mother: " + mother_id + ", father: " + father_id + ". " + str(error))
            
            raise(error)
        
        return variants
    
    def add_single_variant(self, variants, var, gender, line):
        """ adds a single variant to a vcf dictionary indexed by position key
        
        Args:
            variants: list of variants for an individual
            var: single Variant object
            gender: gender of the individual (eg "M"/"F" or "1"/"2")
            line: list of elements of the VCF line for the variant
        """
        
        # Complete the variant setup, now that the variant has passed the
        # filtering. If we do this earlier, it slows all the unneeded variants.
        var.add_format(line[8], line[9])
        var.add_vcf_line(line)
        var.set_gender(gender)
        
        try:
            var.set_genotype()
            variants.append(var)
        except ValueError:
            # we only get ValueError when the genotype cannot be set, which
            # occurs for x chrom male heterozygotes (an impossible genotype)
            if var.get_chrom() == var.debug_chrom and var.get_position() == var.debug_pos:
                print("failed as heterozygous genotype in male on chrX")
            pass
    
    def include_variant(self, line, child_variants, gender, mnvs):
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
            key = (line[0], int(line[1]))
            if key in self.child_keys:
                use_variant = True
            elif line[4] == "<DUP>" or line[4] == "<DEL>":
                var = construct_variant(line, gender, self.known_genes)
                if self.cnv_matcher.has_match(var):
                    use_variant = True
        else:
            var = construct_variant(line, gender, self.known_genes, mnvs)
            if var.passes_filters():
                use_variant = True
            
        return use_variant
        
    def open_individual(self, individual, child_variants=False, mnvs=None):
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
        
        if individual is None:
            return []
        
        path = individual.get_path()
        logging.info("sample path: {}".format(path))
        gender = individual.get_gender()
        
        # open the vcf, and adjust the position in the file to immediately after
        # the header, so we can run through the variants
        vcf = open_vcf(path)
        exclude_header(vcf)
        
        variants = []
        for line in vcf:
            line = line.strip().split("\t")
            
            # check if we want to include the variant or not
            if self.include_variant(line, child_variants, gender, mnvs):
                var = construct_variant(line, gender, self.known_genes, mnvs)
                self.add_single_variant(variants, var, gender, line)
        
        vcf.close()
        
        return variants
    
    def load_trio(self, family):
        """ opens and parses the VCF files for members of the family trio.
        
        We need to load the VCF data for each of the members of the trio. As a
        bare minimum we need VCF data for the child in the family. Occasionally
        we lack parents for the child, so we create blank entries when that
        happens.
        """
        
        mnvs = get_mnv_candidates(family.child.get_path())
        
        # open the childs VCF file, and get the variant keys, to check if they
        # are in the parents VCF
        child = self.open_individual(family.child, mnvs=mnvs)
        self.child_keys = set([var.get_key() for var in child])
        
        self.child_header = get_vcf_header(family.child.get_path())
        self.cnv_matcher = MatchCNVs(child)
        
        mother = self.open_individual(family.mother, child_variants=True)
        father = self.open_individual(family.father, child_variants=True)
        
        return self.combine_trio_variants(family, child, mother, father)
    
    def combine_trio_variants(self, family, child_vars, mother_vars, father_vars):
        """ for each variant, combine the trio's genotypes into TrioGenotypes
        
        Args:
            child_vars: list of Variant objects for the child
            mother_vars: list of Variant objects for the mother
            father_vars: list of Variant objects for the father
        
        Returns:
            list of TrioGenotypes objects for the family
        """
        
        mom_cnvs = MatchCNVs(mother_vars)
        dad_cnvs = MatchCNVs(father_vars)
        
        variants = []
        for child in child_vars:
            
            mother, father = None, None
            if family.has_parents():
                mother = self.get_parental_var(child, mother_vars, family.mother.get_gender(), mom_cnvs)
                father = self.get_parental_var(child, father_vars, family.father.get_gender(), dad_cnvs)
            
            trio = TrioGenotypes(child.get_chrom(), child.get_position(), SNV.debug_chrom, SNV.debug_pos)
            trio.add_child(child)
            trio.add_mother(mother)
            trio.add_father(father)
            
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
        default_ref = '0/0'
        
        # if the variant is a CNV, the corresponding variant might not match
        # the start site, so we look a variant that overlaps
        if var.is_cnv() and matcher.has_match(var):
            key = matcher.get_overlap_key(key)
            
        for parental in parental_vars:
            if key == parental.get_key():
                return parental
        
        # if the childs variant does not exist in the parents VCF, then we
        # create a default variant for the parent
        if var.is_cnv():
            parental = CNV(var.chrom, var.position, var.variant_id, var.ref_allele, '<REF>', var.filter)
            default_ref = 'REF'
            parental.add_info('END=1000000000')
        else:
            parental = SNV(var.chrom, var.position, var.variant_id, var.ref_allele, ','.join(var.alt_alleles), var.filter)
        
        parental.add_format('GT', default_ref)
        parental.set_gender(gender)
        parental.set_genotype()
        
        return parental
    
    def filter_de_novos(self, variants, pp_filter):
        """ filter the de novos variants in the VCF files
        
        Args:
            variants: list of TrioGenotypes objects.
            pp_filter float between 0 and 1, being the threshold for the PP_DNM filter
        
        Returns:
            a list of TrioGenotypes without the de novo variants that failed the
            de novo filter.
        """
        
        return [ x for x in variants if x.passes_de_novo_checks(pp_filter) ]
