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

class PostInheritanceFilter(object):
    """ Post-inheritance variant filter for clinical filtering code.

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
        
        self.variants = self.filter_polyphen(self.variants)
        
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
                log_str = "{}\t{} dropped from excess CNVs in proband".format(self.family.child.get_id(), var)
                logging.info(log_str)
                if var.get_chrom() == self.debug_chrom and var.get_position() == self.debug_pos:
                    print(log_str)
        
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
                    log_str = "{}\t{} dropped from low MAF in non-biallelic" \
                        "variant".format(self.family.child.get_id(), var)
                    logging.info(log_str)
                    if var.get_chrom() == self.debug_chrom and var.get_position() == self.debug_pos:
                        print(log_str)
        
        return passed_vars
    
    def get_polyphen_for_genes(self, var, hgnc):
        """ get the polyphen predictions for a variant for specific gene symbols
        
        Genetic variants can lie within multiple genes, each of which has a
        polyphen prediction (if the site has an appropriate functional
        consequence). This function pulls out the polyphen predictions for
        specific gene symbols, in order to check whether a candidate variant
        has a "benign" polyphen prediction.
        
        Args:
            var: TrioGenotypes object for a variant
            hgnc: list of gene symbols
        
        Returns:
            list of polyphen predictions
        """
        
        # find the HGNC symbol positions in the partner variant that match the
        # HGNC symbols
        genes = var.get_genes()
        pos = [ x for x in range(len(genes)) if genes[x] is not None and genes[x] in hgnc ]
        
        polyphen = []
        if "PolyPhen" in var.child.info:
            # get the polyphen predictions for the genes matching the required
            # gene symbols. NOTE: This does not account for multi-allelic sites,
            # but we don't examine compound hets at multi-allelic sites, so this
            # shouldn't be a problem.
            polyphen = [ var.child.info["PolyPhen"].split("|")[n] for n in pos ]
            
            # remove the numeric scores from the annotations
            polyphen = [ x.split("(")[0] for x in polyphen ]
        
        return polyphen
    
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
        
        for (var, check, inh, hgnc) in variants:
            
            # check if the variant on it's own would pass
            passes = "benign" not in self.get_polyphen_for_genes(var, hgnc) or \
                    var.get_trio_genotype() == var.get_de_novo_genotype()
            
            # check all of the other variants to see if any are in the same
            # gene, compound_het, and polyphen benign
            benign_matches = [ self.has_compound_match(var, x, variants) for x in hgnc ]
            
            # exclude HGNC symbols where partner variants are polyphen benign
            hgnc = [ hgnc[x] for x in range(len(hgnc)) if not(benign_matches[x]) ]
            
            # check if any of the genes for the partner variants have damaging
            # consequences
            benign_match = not(any([ not x for x in benign_matches ]))
            
            if passes and not benign_match:
                passed_vars.append((var, check, inh, hgnc))
            else:
                log_str = "{}\t{} dropped from polyphen prediction".format(self.family.child.get_id(), var)
                logging.info(log_str)
                if var.get_chrom() == self.debug_chrom and var.get_position() == self.debug_pos:
                    print(log_str)
        
        return passed_vars
    
    def has_compound_match(self, var, hgnc, variants):
        """ for a compound var, find if its partner is also polyphen benign
        
        Check all of the other variants to see if any are in the same
        gene, compound_het, and polyphen benign.
        
        Args:
            var: TrioGenotypes object
            hgnc: HGNC symbol that we need to match for the partner variant
            variants: list of (variant, check, inheritance, gene) tuples
        
        Returns:
            True/false for whether there is a compound het match
        """
        
        # get a list of the variants in the gene
        compound_vars = [ x[0] for x in variants if hgnc in x[3] and "compound_het" in x[1] ]
        
        if len(compound_vars) == 0:
            return False
        
        # run through the variants, find all the variants that are not benign,
        # or are benign but de novo.
        not_benign = []
        for alt_var in compound_vars:
            if alt_var.get_trio_genotype() == alt_var.get_de_novo_genotype():
                not_benign.append(alt_var)
            elif "benign" not in self.get_polyphen_for_genes(alt_var, hgnc):
                not_benign.append(alt_var)
        
        # if we have more than two non-benign variants with different genotypes,
        # then we don't want to exclude these variants. Unless we lack parents
        # since those will have the same trio genotypes by virtue of having
        # "NA" values for the parental genotypes.
        genotypes = set([ x.get_trio_genotype() for x in not_benign ])
        if "NA" in var.get_trio_genotype():
            genotypes = [ x.get_trio_genotype() for x in not_benign ]
        
        return len(genotypes) <= 1
    
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
                log_str = "{}\t{} dropped from ExAC frequency count".format(self.family.child.get_id(), var)
                logging.info(log_str)
                if var.get_chrom() == self.debug_chrom and var.get_position() == self.debug_pos:
                    print(log_str)
            
            if inh != []:
                passed_vars.append((var, check, inh, hgnc))
        
        return passed_vars
