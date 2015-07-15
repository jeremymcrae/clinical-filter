""" class for parsing VCF INFO fields
"""

class VariantInfo(object):
    """ parses the VCF info field
    """
    
    # Here are the VEP consequences, ranked in severity from the most severe to
    # the least severe as defined at:
    # http://www.ensembl.org/info/genome/variation/predicted_data.html
    severity = {"transcript_ablation": 0, "splice_donor_variant": 1, \
        "splice_acceptor_variant": 2, "stop_gained": 3, "frameshift_variant": 4, \
        "stop_lost": 5, "start_lost": 6, "initiator_codon_variant": 7, \
        "inframe_insertion": 8, "inframe_deletion": 9, "missense_variant": 10, \
        "protein_altering_variant": 11, "transcript_amplification": 12, \
        "splice_region_variant": 13, "incomplete_terminal_codon_variant": 14, \
        "synonymous_variant": 15, "stop_retained_variant": 16, \
        "coding_sequence_variant": 17, "mature_miRNA_variant": 18, \
        "5_prime_UTR_variant": 19, "3_prime_UTR_variant": 20, \
        "intron_variant": 21, "NMD_transcript_variant": 22, \
        "non_coding_exon_variant": 23, "non_coding_transcript_exon_variant": 24, \
        "nc_transcript_variant": 25, "upstream_gene_variant": 26, \
        "downstream_gene_variant": 27, "TFBS_ablation": 28, \
        "TFBS_amplification": 29, "TF_binding_site_variant": 30,
        "regulatory_region_variant": 31, "regulatory_region_ablation": 32,
        "regulatory_region_amplification": 33, "feature_elongation": 34,
        "feature_truncation": 35, "intergenic_variant": 36}
    
    # define the set of loss-of-function consequences
    lof_consequences = set(["transcript_ablation", "splice_donor_variant", \
        "splice_acceptor_variant", "stop_gained", "frameshift_variant",  \
        "coding_sequence_variant"])
    
    # define the set of missense (or non loss-of-function) consequences
    missense_consequences = set(["stop_lost", "initiator_codon_variant", \
        "inframe_insertion", "inframe_deletion", "missense_variant", \
        "transcript_amplification", "start_lost", "protein_altering_variant"])
    
    # define the populations who have minor allele frequencies in the INFO
    populations = set(["AFR_AF", "AMR_AF", "ASN_AF", "DDD_AF", "EAS_AF", \
        "ESP_AF", "EUR_AF", "MAX_AF", "SAS_AF", "UK10K_cohort_AF"])
    
    # create static variables (which will be set externally before any class
    # objects are created)
    known_genes = None
    debug_chrom = None
    debug_pos = None
    
    def add_info(self, info_values):
        """Parses the INFO column from VCF files.
        
        Args:
            info_values: INFO text from a line in a VCF file
            tags: the tags dict
        """
        
        for item in info_values.split(";"):
            if "=" in item:
                try:
                    key, value = item.split("=")
                except ValueError:
                    pos = item.index("=")
                    key = item[:pos]
                    value = item[pos + 1:]
            else:
                key, value = item, True
            self.info[key] = value
        
        self.set_gene_from_info()
        self.set_consequence()
    
    def has_info(self):
        """ checks if the INFO field has been parsed and added to the object
        """
        
        return self.info != {}
    
    def get_range(self):
        """ gets the range for the CNV
        """
        
        start_position = self.get_position()
        
        if self.is_cnv():
            end_position = start_position + 10000
            if self.has_info() and "END" in self.info:
                end_position = int(self.info["END"])
        else:
            end_position = start_position
        
        return (start_position, end_position)
    
    def set_gene_from_info(self):
        """ sets a gene to the var using the info. CNVs and SNVs act differently
        """
        
        self.gene = None
        
        if "HGNC" in self.info:
            self.gene = []
            for pos in range(len(self.info["HGNC"].split(","))):
                self.gene.append(self.get_genes_for_allele(pos))
            
            # pull out the gene list for single allele variants
            if len(self.gene) == 1:
                self.gene = self.gene[0]
        
        # some genes lack an HGNC entry, but do have an HGNC_ALL entry. The
        # HGNC_ALL entry is a "&"-separated list of Vega symbols.
        elif self.gene is None and "HGNC_ALL" in self.info:
            self.gene = self.info["HGNC_ALL"].split("&")
        # If we are not using a set of known genes, we still want to check
        # variants that haven't been annotated with a HGNC, since some of these
        # have a functional VEP annotation, presumably due to difficulties in
        # identifying an HGNC symbol. We don't need to worry about this when
        # using a set of known genes, since those should all have HGNC symbols.
        elif self.gene is None and self.known_genes is None:
            self.gene = "{0}:{1}".format(self.chrom, self.position)
    
    def get_genes_for_allele(self, position):
        """ gets list of gene symbols for an allele, prioritising HGNC symbols.
        
        We have a variety of gene symbol sources for a variant. The INFO field
        can contain HGNC, SYMBOL, ENSG, ENST, ENSP and ENSR. HGNC is HGNC gene
        symbol, SYMBOL is VEGA-derived gene symbol, ENSG is Ensembl gene ID,
        ENST is Ensembl transcript ID, ENSP is Ensembl protein ID, and ENSR is
        Ensembl regulatory ID.
        
        In order for variants to check for compound heterozygotes, we match
        variants by gene symbols. Ideally each variant would have a list of HGNC
        symbols that it occurs in, but this is not the case. Many variants have
        HGNC entries such as ".|GENE1|GENE2|.", where the "." indicates no
        symbol available. We fill in the missing symbols by successively looking
        through the SYMBOL field (from VEGA-derived symbols), then the ENSG IDs
        and so on. Some entries will always lack a value, such as for variants
        in transcription factor binding sites, where no symbol will ever be
        available for the entry.
        
        Args:
            position: integer position for the allele to be examined within a
                comma-separated list.
        
        Returns:
            list of gene symbols, filled in from the HGNC, SYMBOL, ENSG fields
            where available.
        """
        
        # set the list of fields to check, in order of their priority.
        fields = ["HGNC", "SYMBOL", "ENSG", "ENST", "ENSP", "ENSR"]
        fields = [ x for x in fields if x in self.info ]
        
        genes = None
        for field in fields:
            symbols = self.info[field].split(",")[position].split("|")
            if genes is None:
                genes = symbols
            
            # find which positions of the gene list are missing a symbol
            blanks = [ x for x,item in enumerate(genes) if item in ["", "."] ]
            
            # for the missing symbol positions, if the current list of symbols
            # contains a value, swap that value into the gene list.
            for x in blanks:
                if symbols[x] != "":
                    genes[x] = symbols[x]
        
        if genes is not None:
            genes = [ x if x not in ["", "."] else None for x in genes ]
        
        return genes
    
    def get_genes(self):
        """ split a gene string into list of gene names
        
        Returns:
            list of gene IDs
        """
        
        genes = self.gene
        if self.gene is None:
            genes = []
        
        return genes
    
    def set_consequence(self):
        """ makes sure a consequence field is available in the info dict
        """
        
        cq = None
        if "CQ" in self.info:
            cq = self.info["CQ"].split("|")
        
        if "," in self.alt_allele:
            (cq, self.gene, enst) = self.correct_multiple_alt(cq)
            
            if "ENST" in self.info:
                self.info["ENST"] = enst
        
        self.consequence = cq
    
    def get_per_gene_consequence(self, hgnc_symbol):
        """ find the VEP consequences for a HGNC symbol.
        
        Some variants have consequences annotated for mutlple genes (because
        the variant lies in multiple genes). We occasionally want to know what
        the consequences are for a specific gene. It's possible that a variant
        might have two symbols of "", in which case we return the consequences
        for both instances of this symbol.
        
        Args:
            hgnc_symbol: HGNC symbol for which we wish to check VEP consequence.
        
        Returns:
            List of VEP consequences for the variant. If hgnc_symbol is None,
            then we return all the consequences listed for the variant.
            Otherwise we find the consequences listed for the gene symbol.
        """
        
        if hgnc_symbol is None:
            return self.consequence
        
        # find the positions of the genes that match the curent HGNC symbol.
        # There could be multiple matches if the symbol is "".
        pos = [ i for i, item in enumerate(self.get_genes()) if item == hgnc_symbol ]
        
        # At one point, the VCFs lacked per gene consequences, but could have
        # multiple gene symbols (if they lacked a HGNC field but did have a
        # HGNC_ALL field). These variants will have multiple genes, but only one
        # consequence. Return the consequence as is, in order to retain the
        # same output for those VCFs.
        if len(self.get_genes()) > 1 and len(self.consequence) == 1:
            return self.consequence
        
        return [ self.consequence[n] for n in pos ]
        
    def correct_multiple_alt(self, cq):
        """ gets correct consequence, HGNC and ensembl IDs for multiple alt vars
        
        Some variants have multiple alts, so we need to select the alt with
        the most severe consequence. However, in at least one version of the
        VCFs, one of the alts could have zero depth, which I believe resulted
        from the population based multi-sample calling. We need to drop the
        consequences recorded for zero-depth alternate alleles before finding
        the most severe.
        
        We also correct HGNC and Ensembl transcript IDs, since otherwise the
        multiple alt variants will contain as many repeated instances of the
        symbols as there are alt alleles.
        
        Args:
            cq: list of VEP annotated consequences
        
        Returns:
            tuple of ([consequence], [HGNC symbol], and [ensembl transcript ID])
            lists
        """
        
        # we join the consequence list (which has been split by gene), so that
        # splitting by allele can take priority
        cq = "|".join(cq)
        
        # get the consequence string, then split the entries for each allele by
        # gene
        cq = cq.split(",")
        cq = [x.split("|") for x in cq]
        
        # drop the consequence for zero-depth alleles
        if "AC" in self.info:
            # check if alleles are present in the individual
            present = [ x != "0" for x in self.info["AC"].split(",") ]
            
            # exclude the consequences for alleles with zero alleles
            cq[:] = [ item for i,item in enumerate(cq) if present[i] ]
    
        cq = self.get_most_severe_consequence(cq)
        
        hgnc = self.get_genes()
        if len(hgnc) > 0:
            lengths = [ len(x) for x in hgnc ]
            max_len = max(lengths)
            pos = lengths.index(max_len)
            hgnc = hgnc[pos]
        
        enst = None
        if "ENST" in self.info:
            enst = self.info["ENST"].split(",")
            enst = [x.split("|") for x in enst]
            enst = "|".join(enst[0])
        
        return (cq, hgnc, enst)
    
    def get_most_severe_consequence(self, consequences):
        """ get the most severe consequence from a list of vep consequence terms
        
        Args:
            consequences: list of VEP consequence strings, or list of lists
        
        Returns:
            the most severe consequence string
        """
        
        # If we have passed in a list of lists, such as for multiple alleles,
        # then we consolidate all of the consequences into per gene lists
        # (rather than per allele lists), and check for the most severe
        # consequence within each geen list.
        if type(consequences[0]) is list:
            new_list = []
            max_len = max([ len(x) for x in consequences ])
            for i in range(max_len):
                per_gene = [ x[i] for x in consequences if i < len(x)]
                new_list.append(self.get_most_severe_consequence(per_gene))
            
            return new_list
        
        most_severe = ""
        most_severe_score = 1000
        for cq in consequences:
            if self.severity[cq] < most_severe_score:
                most_severe = cq
                most_severe_score = self.severity[cq]
        
        return most_severe
       
    def is_lof(self, hgnc_symbol=None):
        """ checks if a variant has a loss-of-function consequence
        
        Args:
            hgnc_symbol: HGNC symbol for which we wish to check VEP consequence.
                By default we check all the consequences listed for the variant.
        """
        
        if self.consequence is None:
            return False
        
        cq = self.get_per_gene_consequence(hgnc_symbol)
        
        return len(set(cq) & self.lof_consequences) > 0
    
    def is_missense(self, hgnc_symbol=None):
        """ checks if a variant has a missense-styled consequence
        
        Args:
            hgnc_symbol: HGNC symbol for which we wish to check VEP consequence.
                By default we check all the consequences listed for the variant.
        """
        
        if self.consequence is None:
            return False
        
        cq = self.get_per_gene_consequence(hgnc_symbol)
        
        return len(set(cq) & self.missense_consequences) > 0
    
    def get_allele_frequency(self, values):
        """ extracts the allele frequency float from a VCF string
        
        The allele frequency for a population can be encoded in several ways,
        either as a single float (eg "0.01"), or as a missing value (eg "."),
        or there can be a list of allele frequencies for the different alternate
        alleles for the variant (eg "0.01,0.05,0.06"), or list containing floats
        and missing values. We need to return the allele frequency as a float,
        but if there are multiple allele frequencies, we return the largest
        float.
        
        Args:
            values: string for allele frequency eg "0.01" or ".", or
                "0.01,.,0.06". Sometimes we might even get values passed in as
                a float, or a None type.
        
        Returns:
            allele frequency as float, or None, if no frequency available
        """
        
        if isinstance(values, float):
            return values
        
        if values is None:
            return None
        
        values = values.split(",")
        values = [ float(x) for x in values if self.is_number(x) ]
        
        if values == []:
            return None
        
        return max(values)
    
    def is_number(self, value):
        """ determines whether a value represents a number.
        
        Sometimes the MAF reported for a variant is ".", or even ".,.", which
        are not numbers and are in fact NA values, but would cause the variant
        not to pass the MAF filter. instead check if the value can be
        converted to a float.
        
        Args:
            value: a string or other number
        
        Returns:
            True or False for whether the value can be converted to a float.
        """
        
        if value is None:
            return False
        
        try:
            value = float(value)
            return True
        except ValueError:
            return False
        
        return False
    
    def find_max_allele_frequency(self):
        """gets the maximum allele frequency for a variant in a VCF record
        
        Finds the maximum allele frequency recorded for a variant across
        different populations.
        
        Args:
            populations: list of population IDs to search
          
        Returns:
            the maximum allele frequency found within the populations in the
            variant record
        """
        
        max_freq = None
        # check all the populations with MAF values recorded for the variant
        # (typically the 1000 Genomes populations (AFR_AF, EUR_AF etc), any
        # internal population (e.g. DDD_AF), and a MAX_AF field)
        for key in set(self.populations) & set(self.info):
            frequency = self.get_allele_frequency(self.info[key])
            if frequency is None:
                continue
            
            if max_freq is None or frequency > max_freq:
                max_freq = frequency
        
        return max_freq
    
    
