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

import logging

class PostInheritanceFilter(object):
    """ Apply some post inheritance filters to flagged variants
    """
    
    def __init__(self, variants, family, debug_chrom=None, debug_pos=None):
        """intialise the class with the some definitions
        """
        
        self.variants = variants
        self.debug_chrom = debug_chrom
        self.debug_pos = debug_pos
        self.family = family
    
    def filter_variants(self):
        """ loads trio variants, and screens for candidate variants
        """
        
        # if we have flagged CNVs on three different chroms, drop all CNVs,
        # since the sample is sufficiently anomalous
        if self.count_cnv_chroms(self.variants) > 2:
            self.variants = self.remove_cnvs(self.variants)
        
        # and filter by a lower MAF threshold
        self.variants = self.filter_by_maf(self.variants)
        
        self.variants = self.filter_exac(self.variants)
        
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
        for (var, check, inh, hgnc) in variants:
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
        for (var, check, inh, hgnc) in variants:
            if not var.is_cnv():
                passed_vars.append((var, check, inh, hgnc))
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
        for (var, check, inh, hgnc) in variants:
            max_maf = var.child.find_max_allele_frequency()
            if max_maf is None: # set maf=NA to 0 to reduce later checks
                max_maf = 0
            
            if inh == ["Biallelic"]:
                passed_vars.append((var, check, inh, hgnc))
            # variants with multiple inheritance types should be left as
            # Biallelic if the other inheritance type fails the MAF threshold
            elif "Biallelic" in inh and max_maf >= 0.001:
                passed_vars.append((var, check, ["Biallelic"], hgnc))
            else:
                if max_maf <= 0.001 and self.family.has_parents():
                    passed_vars.append((var, check, inh, hgnc))
                elif max_maf <= 0.0001 and not self.family.has_parents():
                    passed_vars.append((var, check, inh, hgnc))
                else:
                    logging.debug(str(var) + " dropped from low MAF in non-biallelic variant")
                    if var.get_chrom() == self.debug_chrom and var.get_position() == self.debug_pos:
                        print(str(var) + " dropped from low MAF in non-biallelic variant")
        
        return passed_vars
    
    def filter_exac(self, variants):
        """ drop variants based on ExAC frequencies
        
        The frequency threshold depends on the inheritance mode.
        
        Args:
            variants: list of (variant, check, inheritance, gene) tuples
        
        Returns:
            returns list of tuples without variants where the variant is on
            has too high a frequency in ExAC.
        """
        
        passed_vars = []
        
        for (var, check, inh, hgnc) in variants:
            
            # figure out what the het and hemi counts are in ExAC (if available)
            hemi, het = 0, 0
            if "AC_Hemi" in var.child.info and var.get_chrom() == "X":
                hemi = sum([int(x) for x in var.child.info["AC_Hemi"].split(",")])
            if "AC_Het" in var.child.info:
                het = sum([int(x) for x in var.child.info["AC_Het"].split(",")])
            
            geno = var.get_trio_genotype()
            # filter out hemizygous variants on chrX in males. Autosomal
            # and female chrX variants should pass through unfiltered.
            # We don't filter out de novo variants based on the ExAC hemizygous
            # count. We only apply this filter to inherited variants.
            if "Hemizygous" in inh and self.family.child.is_male() and hemi > 0 and \
                geno != var.get_de_novo_genotype() and geno[1:] != ("NA", "NA"):
                    inh.remove("Hemizygous")
            
            # filter out monoallelic variants with high ExAC het counts.
            if var.get_chrom() != "X" and het > 4 and "Monoallelic" in inh:
                inh.remove("Monoallelic")
            elif var.get_chrom() == "X" and (het + hemi) > 4 and "X-linked dominant" in inh:
                inh.remove("X-linked dominant")
            
            if inh == []:
                logging.debug("{} dropped from ExAC frequency count".format(var))
                if var.get_chrom() == self.debug_chrom and var.get_position() == self.debug_pos:
                    print("{} dropped from ExAC frequency count".format(var))
            
            if inh != []:
                passed_vars.append((var, check, inh, hgnc))
        
        return passed_vars
