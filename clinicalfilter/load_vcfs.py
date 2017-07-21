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

from clinicalfilter.variant.info import Info
from clinicalfilter.variant.variant import Variant
from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.trio_genotypes import TrioGenotypes
from clinicalfilter.utils import open_vcf, get_vcf_header, exclude_header, \
    construct_variant
from clinicalfilter.multinucleotide_variants import get_mnv_candidates

class LoadVCFs(object):
    """ load VCF files for a trio
    """
    
    def __init__(self, probands_n, maf_tags, known_genes, last_base, debug_chrom, debug_pos):
        """ intitalise the class with the filters and tags details etc
        
        Args:
            probands_n: count of how many probands are to be analysed
            maf_tags: list of populations who have minor allele frequencies in
                the INFO.
            known_genes: dictionary of genes known to be involved with genetic
                disorders.
            last_base: set of sites in genome at conserved last base of exons,
                where we upgrade the severity of variants to loss-of-function.
            debug_chrom: chromosome string, to give more information about why
                a variant fails to pass the filters.
            debug_pos: chromosome position, to give more information about why
                a variant fails to pass the filters.
        """
        
        self.counter = 0
        self.probands_n = probands_n
        
        # define several parameters of the variant classes, before we have
        # initialised any class objects
        SNV.set_known_genes(known_genes)
        CNV.set_known_genes(known_genes)
        Info.set_last_base_sites(last_base)
        Info.set_populations(maf_tags)
        
        SNV.set_debug(debug_chrom, debug_pos)
        CNV.set_debug(debug_chrom, debug_pos)
    
    def get_trio_variants(self, family, pp_filter):
        """ loads the variants for a trio
        
        Args:
            family: Family object for a trio
            pp_filter float between 0 and 1, being the threshold for the PP_DNM filter
        
        Returns:
            list of filtered variants for a trio, as TrioGenotypes objects
        """
        
        self.counter += 1
        logging.info("opening trio {} of {}".format(self.counter, self.probands_n))
        
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
    
    def include_variant(self, line, child_variants, gender, mnvs):
        """ check if we want to include the variant or not
        
        Args:
            line: list of elements from the VCF line for the variant.
            child_variants: list of variants that passed in the child, so we can
                quickly assess parental variants. This is None when screening
                the child.
            gender: the gender of the proband (used in CNV filtering).
            mnvs: dictionary of (chrom, pos), MNV_code pairs for known
                multinucleotide variant sites  within the proband.
        
        Returns:
            True/False for whether to include the variant.
        """
        
        if child_variants is not None:
            key = (line[0], int(line[1]))
            return key in child_variants
        
        var = construct_variant(line, gender, mnvs)
        return var.passes_filters()
        
    def open_individual(self, individual, child_variants=None, mnvs=None):
        """ Convert VCF to TSV format. Use for single sample VCF file.
        
        Obtains the VCF data for a single sample. This function optionally
        filters the lines of the VCF file that pass defined criteria, in order
        to reduce memory usage.
        
        Args:
            individual: Person object for individual
            child_variants: True/False for whether variants have been filtered
                for the proband (if so, we can simply check the parent's
                variants for matches in the child's variants).
            mnvs: dictionary
        
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
            
            try:
                # check if we want to include the variant or not
                if self.include_variant(line, child_variants, gender, mnvs):
                    var = construct_variant(line, gender, mnvs)
                    var.add_vcf_line(line)
                    variants.append(var)
            except ValueError:
                # we only get ValueError when the genotype cannot be set, which
                # occurs for x chrom male heterozygotes (an impossible genotype)
                if line[0] == SNV.debug_chrom and int(line[1]) == SNV.debug_pos:
                    print("failed as heterozygous genotype in male on chrX")
                continue
        
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
        keys = set([var.get_key() for var in child])
        
        self.child_header = get_vcf_header(family.child.get_path())
        
        mother = self.open_individual(family.mother, child_variants=keys)
        father = self.open_individual(family.father, child_variants=keys)
        
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
        
        variants = []
        for child in child_vars:
            
            mom, dad = None, None
            if family.has_parents():
                mom = self.get_parental_var(child, mother_vars, family.mother)
                dad = self.get_parental_var(child, father_vars, family.father)
            
            trio = TrioGenotypes(child.get_chrom(), child.get_position(),
                child, mom, dad, SNV.debug_chrom, SNV.debug_pos)
            
            variants.append(trio)
        
        return variants
    
    def get_parental_var(self, var, parental_vars, parent):
        """ get the corresponding parental variant to a childs variant, or
        create a default variant with reference genotype.
        
        Args:
            var: childs var, as Variant object
            parental_vars: list of parental variants
            parent: Person object for the parent
        
        Returns:
            returns a Variant object, matched to the proband's variant
        """
        
        key = var.get_key()
        
        for parental in parental_vars:
            if not var.is_cnv() and key == parental.get_key():
                return parental
        
        # if the childs variant does not exist in the parents VCF, then we
        # create a default variant for the parent
        Var = SNV
        keys, sample = 'GT', '0/0'
        alts = ','.join(var.alt_alleles)
        
        if var.is_cnv():
            Var = CNV
            inh = var.get_cnv_inheritance()
            alts = ("<REF>", )
            if parent.is_male() and inh in ['paternal', 'biparental']:
                alts = var.alt_alleles
            elif parent.is_female() and inh in ['maternal', 'biparental']:
                alts = var.alt_alleles
            
            alts = ','.join(alts)
            # we need to set a format value, so CNV genotypes get set correctly
            keys, sample = 'INHERITANCE', 'uncertain'
        
        return Var(var.chrom, var.position, var.variant_id, var.ref_allele,
            alts, var.qual, var.filter, str(var.info), keys, sample, 
            parent.get_gender())
    
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
