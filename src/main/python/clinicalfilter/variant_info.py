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
        "stop_lost": 5, "initiator_codon_variant": 6, "inframe_insertion": 7, \
        "inframe_deletion": 8, "missense_variant": 9, \
        "transcript_amplification": 10, "splice_region_variant": 11, \
        "incomplete_terminal_codon_variant": 12, "synonymous_variant": 13, \
        "stop_retained_variant": 14, "coding_sequence_variant": 15, \
        "mature_miRNA_variant": 16, "5_prime_UTR_variant": 17, \
        "3_prime_UTR_variant": 18, "intron_variant": 19, \
        "NMD_transcript_variant": 20, "non_coding_exon_variant": 21, \
        "nc_transcript_variant": 22, "upstream_gene_variant": 23, \
        "downstream_gene_variant": 24, "TFBS_ablation": 25, \
        "TFBS_amplification": 26, "TF_binding_site_variant": 27, \
        "regulatory_region_variant": 28, "regulatory_region_ablation": 29, \
        "regulatory_region_amplification": 30, "feature_elongation": 31, \
        "feature_truncation": 32, "intergenic_variant": 33}
    
    # define the set of loss-of-function consequences
    lof_consequences = set(["transcript_ablation", "splice_donor_variant", \
        "splice_acceptor_variant", "stop_gained", "frameshift_variant",  \
        "coding_sequence_variant"])
    
    # define the set of missense (or non loss-of-function) consequences
    missense_consequences = set(["stop_lost", "initiator_codon_variant", \
        "inframe_insertion", "inframe_deletion", "missense_variant", \
        "transcript_amplification"])
    
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
        
        self.set_consequence()
        self.set_gene_from_info()
    
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
        # sometimes the variant lacks an HGNC field
        if "HGNC" in self.info:
            self.gene = self.info["HGNC"]
        
    def set_gene_from_known_gene_overlap(self):
        """ sets the gene according to overlap with the positions of known genes
        """
        
        # don't change anything if we don't have a set of known genes
        if self.known_genes is None:
            return
        
        overlapping = self.get_overlapping_known_genes()
        
        previous = []
        if self.gene is not None:
            previous = self.gene.split(",")
        
        self.gene = ",".join(sorted(set(overlapping + previous)))
    
    def get_overlapping_known_genes(self):
        """ finds the names of known genes that a variant overlaps
        """
        
        (start, end) = self.get_range()
        
        if self.known_genes is None:
            raise ValueError("we don't have a set of known genes to look through")
        
        current_chrom = self.get_chrom()
        current_pos = self.get_position()
        
        overlapping = []
        for gene in self.known_genes:
            if self.known_genes[gene]["chrom"] != current_chrom:
                continue
            
            gene_start = self.known_genes[gene]["start"]
            gene_end = self.known_genes[gene]["end"]
            
            if start <= gene_end and end >= gene_start:
                overlapping.append(gene)
        
        return overlapping
    
    def set_consequence(self):
        """ makes sure a consequence field is available in the info dict
        """
        
        if "CQ" not in self.info:
            self.info["CQ"] = None
        
        cq = self.info["CQ"]
        
        if "," in self.alt_allele:
            (cq, hgnc, enst) = self.correct_multiple_alt(cq)
            
            if "HGNC" in self.info:   
                self.info["HGNC"] = hgnc
                self.info["ENST"] = enst
        
        self.consequence = cq
        
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
            cq: string of VEP annotated comma-separated consequences
        
        Returns:
            tuple of (consequence, HGNC symbol, and ensembl transcript ID)
        """
        
        # get the consequence string
        cq = cq.split(",")
        
        enst = None
        hgnc = None
        if "HGNC" in self.info:
            enst = self.info["ENST"].split(",")
            hgnc = self.info["HGNC"].split(",")
        
        # drop the consequence for zero-depth alleles
        if "AC" in self.info:
            # check if alleles are present in the individual
            present = [ x != "0" for x in self.info["AC"].split(",") ]
            
            # exclude the consequences for alleles with zero alleles
            cq[:] = [ item for i,item in enumerate(cq) if present[i] ]
            
            if "HGNC" in self.info:
                hgnc[:] = [ item for i,item in enumerate(hgnc) if present[i] ]
                enst[:] = [ item for i,item in enumerate(enst) if present[i] ]
        
        cq = self.get_most_severe_consequence(cq)
        
        if "HGNC" in self.info:   
            hgnc = ",".join(sorted(set(hgnc)))
            enst = ",".join(sorted(set(enst)))
        
        return (cq, hgnc, enst)
    
    def get_most_severe_consequence(self, consequences):
        """ get the most severe consequence from a list of vep consequence terms
        
        Args:
            consequences: list of VEP consequence strings
        
        Returns:
            the most severe consequence string
        """
        
        most_severe = ""
        most_severe_score = 1000
        for cq in consequences:
            if self.severity[cq] < most_severe_score:
                most_severe = cq
                most_severe_score = self.severity[cq]
        
        return most_severe
       
    def is_lof(self):
        """ checks if a variant has a loss-of-function consequence
        """
        
        return self.consequence in self.lof_consequences
    
    def is_missense(self):
        """ checks if a variant has a missense-styled consequence
        """
        
        return self.consequence in self.missense_consequences
    
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
    
    
