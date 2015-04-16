""" Postinheritance variant filter for clinical filtering code.

Most variants can be excluded on the basis of their individual characteristics
(minor allele frequency, quality scores etc), but some variants need to be
assessed after the standard filters
    we fail CNVs where an individual has CNVs that pass the filters on three or
        more chromosomes (so we can only filter these out once we have found the
        CNVs that pass the filters).
    We fail CNVs with minor allele frequencies greater than 0.1%, EXCEPT for
        biallelic variants, they have a 1% threshold. We can only filter these
        out once we have found variants that pass the different inheritance
        models - then we can check if they are biallelelic or not.
    We fail SNVs with polyphen=benign, even if a compound pair of the SNV is
        polyphen=benign, we require both SNVs to be polyphen=not benign
"""

import sys
import logging

class PostInheritanceFilter(object):
    """ Apply some post inheritance filters to flagged variants
    """
    
    def __init__(self, variants, debug_chrom=None, debug_pos=None):
        """intialise the class with the some definitions
        """
        
        self.variants = variants
        self.debug_chrom = debug_chrom
        self.debug_pos = debug_pos
    
    def filter_variants(self):
        """ loads trio variants, and screens for candidate variants
        """
        
        # if we have flagged CNVs on three different chroms, drop all CNVs,
        # since the sample is sufficiently anomalous
        if self.count_cnv_chroms(self.variants) > 2:
            self.variants = self.remove_cnvs(self.variants)
        
        # and filter by a lower MAF threshold
        self.variants = self.filter_by_maf(self.variants)
        
        self.variants = self.filter_polyphen(self.variants)
        
        self.filter_exac_hemizygous(self.variants)
        
        return self.variants
    
    def count_cnv_chroms(self, variants):
        """ count the number of different chroms that CNVs are on
        
        Args:
            variants: list of (variant, check, inheritance) tuples
        
        Returns:
            integer count of distinct chromosomes containing CNVs
        """
        
        # count the flagged CNVs
        cnv_chroms = set([])
        for (var, check, inh) in variants:
            if var.is_cnv():
                cnv_chroms.add(var.get_chrom())
        
        return len(cnv_chroms)
    
    def remove_cnvs(self, variants):
        """ remove CNVs from individuals with too many flagged CNVs
        
        Args:
            variants: list of (variant, check, inheritance) tuples
        
        Returns:
            returns list of tuples without CNV variants
        """
        
        # remove CNVs from the list of flagged variants
        passed_vars = []
        for (var, check, inh) in variants:
            if not var.is_cnv():
                passed_vars.append((var, check, inh))
            else:
                logging.debug(str(var) + " dropped from excess CNVs in proband")
                if var.get_chrom() == self.debug_chrom and var.get_position() == self.debug_pos:
                    print(str(var) + " dropped from excess CNVs in proband")
        
        return  passed_vars
    
    def filter_by_maf(self, variants):
        """ filter for low MAF threshold, except for Biallelic variants
        
        Args:
            variants: list of (variant, check, inheritance) tuples
        
        Returns:
            returns list of tuples without high maf variants
        """
        
        passed_vars = []
        for (var, check, inh) in variants:
            # populations = var.child.tags["MAX_MAF"]
            max_maf = var.child.find_max_allele_frequency()
            if max_maf is None: # set maf=NA to 0 to reduce later checks
                max_maf = 0
            
            if inh == "Biallelic":
                passed_vars.append((var, check, inh))
            # variants with multiple inheritance types should be left as
            # Biallelic if the other inheritance type fails the MAF threshold
            elif "Biallelic" in inh and max_maf >= 0.001:
                passed_vars.append((var, check, "Biallelic"))
            else:
                if max_maf <= 0.001:
                    passed_vars.append((var, check, inh))
                else:
                    logging.debug(str(var) + " dropped from low MAF in non-biallelic variant")
                    if var.get_chrom() == self.debug_chrom and var.get_position() == self.debug_pos:
                        print(str(var) + " dropped from low MAF in non-biallelic variant")
        
        return passed_vars
    
    def filter_polyphen(self, variants):
        """ filter variants based on polyphen predictions
        
        filter out compound hets where both have benign predictions from
        polyphen, but retain compound hets where only one is polyphen benign.
        Also filter out single variants where polyphen predicts benign.
        
        Args:
            variants: list of (variant, check, inheritance) tuples
        
        Returns:
            returns list of tuples without polyphen benign variants
        """
        
        passed_vars = []
        
        for (var, check, inh) in variants:
            passes = False
            
            # check if the variant on it's own would pass
            if "PolyPhen" not in var.child.info or \
                    not var.child.info["PolyPhen"].startswith("benign") or \
                    var.get_trio_genotype() == var.get_de_novo_genotype() or \
                    var.get_trio_genotype()[1:] == ("NA", "NA"):
                passes = True
            
            # check all of the other variants to see if any are in the same
            # gene, compound_het, and polyphen benign
            benign_match = self.has_compound_match(var, variants)
            
            if passes and not benign_match:
                passed_vars.append((var, check, inh))
            else:
                logging.debug(str(var) + " dropped from polyphen prediction")
                if var.get_chrom() == self.debug_chrom and var.get_position() == self.debug_pos:
                    print(str(var) + " dropped from polyphen prediction")
        
        return passed_vars
    
    def has_compound_match(self, var, variants):
        """ for a compound var, find if its partner is also polyphen benign
        
        Check all of the other variants to see if any are in the same
        gene, compound_het, and polyphen benign.
        
        Args:
            var: TrioGenotypes object
            variants: list of (variant, check, inheritance) tuples
        
        Returns:
            True/false for whether there is a compound het match
        """
        
        # get a list of the variants in the gene
        compound_vars = []
        for (alt_var, alt_check, alt_inh) in variants:
            # ignore if we are looking at another gene
            if var.child.gene != alt_var.child.gene:
                continue
            if "compound_het" not in alt_check:
                continue
            compound_vars.append(alt_var)
        
        if len(compound_vars) == 0:
            return False
        
        # run through the variants, find all the variants that are not benign,
        # or are benign but de novo.
        not_benign = []
        for alt_var in compound_vars:
            if "PolyPhen" not in alt_var.child.info:
                not_benign.append(alt_var)
            elif alt_var.get_trio_genotype() == alt_var.get_de_novo_genotype():
                not_benign.append(alt_var)
            elif alt_var.get_trio_genotype()[1:] == ("NA", "NA"):
                not_benign.append(alt_var)
            elif not alt_var.child.info["PolyPhen"].startswith("benign"):
                not_benign.append(alt_var)
        
        # if we have more than two non-benign variants with different genotypes,
        # then we don't want to exclude these variants. Unless we lack parents
        # since those will have the same trio genotypes by virtue of having
        # "NA" values for the parental genotypes.
        genotypes = set([ x.get_trio_genotype() for x in not_benign ])
        if var.get_trio_genotype()[1:] == ("NA", "NA"):
            genotypes = [ x.get_trio_genotype() for x in not_benign ]
        
        return len(genotypes) <= 1
    
    def filter_exac_hemizygous(self, variants):
        """ drop inherited chrX male variants with ExAC hemizygous frequencies >0
        
        Args:
            variants: list of (variant, check, inheritance) tuples
        
        Returns:
            returns list of tuples without variants where the variant is on
            chrX, in a male, inherited, and has a non-zero ExAC hemizygous
            frequency.
        """
        
        passed_vars = []
        
        for (var, check, inh) in variants:
            passes = True
            
            # we only apply this filter to variants on chrX in males. Autosomal
            # and female chrX variants should pass through unfiltered.
            if var.inheritance_type == "XChrMale" and \
                "AC_Hemi" in var.child.info and int(var.child.info["AC_Hemi"]) > 0:
                    passes = False
            
            # we don't filter out de novo variants based on the ExAC hemizygous
            # count. We only apply this filter to inherited variants.
            if var.get_trio_genotype() == var.get_de_novo_genotype():
                passes = True
            elif var.get_trio_genotype()[1:] == ("NA", "NA"):
                passes = True
            
            if passes:
                passed_vars.append((var, check, inh))
        
        return passed_vars
