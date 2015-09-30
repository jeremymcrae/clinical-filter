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

class Info(object):
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
        "start_lost", "initiator_codon_variant",
        "conserved_exon_terminus_variant"])
    
    # define the set of missense (or non loss-of-function) consequences
    missense_consequences = set(["stop_lost", \
        "inframe_insertion", "inframe_deletion", "missense_variant", \
        "transcript_amplification", "protein_altering_variant"])
    
    synonymous_consequences = set(["synonymous_variant"])
    
    # create static variables (set before creating any class instances)
    known_genes = None
    last_base = set([])
    populations = []
    
    @classmethod
    def set_known_genes(cls_obj, known_genes):
        cls_obj.known_genes = known_genes
    
    @classmethod
    def set_last_base_sites(cls_obj, sites):
        cls_obj.last_base = sites
    
    @classmethod
    def set_populations(cls_obj, populations):
        '''define the populations who have minor allele frequencies in the INFO
        '''
        if populations is not None:
            assert type(populations) == list
            cls_obj.populations = populations
    
    def add_info(self, info_values):
        """Parses the INFO column from VCF files.
        
        Args:
            info_values: INFO text from a line in a VCF file
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
        
        masked = self.get_zero_count_alleles(self.info, self.alt_alleles)
        self.genes = self.get_gene_from_info(self.info, self.alt_alleles, masked)
        self.consequence = self.get_consequences(self.info, self.alt_alleles, masked)
    
    def get_info_as_string(self):
        ''' reprocess the info dictionary back into a string, correctly sorted
        '''
        
        info = None
        if self.info != {}:
            info = []
            for key, value in sorted(self.info.items()):
                entry = key
                if value != True:
                    entry = '{}={}'.format(key, value)
                info.append(entry)
            info = ';'.join(info)
        
        return info
    
    def add_info_field(self, key, value):
        ''' add another entry to the info dictionary
        
        Args:
            key: string for the info key
            value: value for the field, or True if the key will be used as a flag
        
        Raises:
            ValueError if the key is already present in the dictionary
        '''
        
        if key not in self.info:
            self.info[key] = value
        else:
            raise ValueError('tried to add an already existing field to the INFO')
    
    def get_range(self):
        """ gets the range for the CNV
        """
        
        start_position = self.get_position()
        
        if self.is_cnv():
            end_position = start_position + 10000
            if "END" in self.info:
                end_position = int(self.info["END"])
        else:
            end_position = start_position
        
        return (start_position, end_position)
    
    def get_gene_from_info(self, info, alts, masked):
        """ sets a gene to the var using the info. CNVs and SNVs act differently
        
        Args:
            info: dictionary of keys and values for the info fields. Contains
                entries for HGNC, SYMBOL, ENSG, ENST, ENSP, ENSR. These are a
                comma-separated list of gene (or transcript, protein etc)
                symbols for the possible alt alleles. Each entry in the
                comma-separated is a pipe-separated list of the genes (or
                transcripts etc) that the particular allele occurs in. CNV (and
                some SNVs in some DDD VCF versions) can also have a HGNC_ALL entry.
                CNVs can also have NUMBERGENES, indicating the number of genes
                that the CNV overlaps, rather than providing the full list of
                HGNC symbols affected.
            alts: list of alternative alleles for the variant
            masked: list of alternative alleles that we don't consider. These
                are identified as alt alleles with zero depth in the individual.
                This can occur due to multi-sample calling.
        
        Returns:
            list of gene lists, one per alternative allele (after removing the
            masked alt alleles.)
        """
        
        genes = None
        
        if "HGNC" in info or "SYMBOL" in info:
            pos = [ i for i, x in enumerate(alts) if x not in masked ]
            genes = [ self.get_genes_for_allele(info, i) for i in pos ]
        # some genes lack an HGNC entry, but do have an HGNC_ALL entry. The
        # HGNC_ALL entry is a "&"-separated list of Vega symbols.
        elif genes is None and "HGNC_ALL" in info:
            genes = [info["HGNC_ALL"].split("&") * (len(alts) - len(masked))]
        elif self.is_cnv() and "HGNC" not in info and "NUMBERGENES" in info:
            if int(info["NUMBERGENES"]) > 0:
                genes = [["."]]
        # If we are not using a set of known genes, we still want to check
        # variants that haven't been annotated with a HGNC, since some of these
        # have a functional VEP annotation, presumably due to difficulties in
        # identifying an HGNC symbol. We don't need to worry about this when
        # using a set of known genes, since those should all have HGNC symbols.
        elif genes is None and self.known_genes is None:
            symbol = '{0}:{1}'.format(self.get_chrom(), self.get_position())
            genes = [[symbol] * (len(alts) - len(masked))]
        
        return genes
    
    def get_genes_for_allele(self, info, position):
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
        fields = [ x for x in fields if x in info ]
        
        genes = None
        for field in fields:
            symbols = info[field].split(",")[position].split("|")
            if genes is None:
                genes = symbols
            
            # find which positions of the gene list are missing a symbol
            blanks = [ x for x,item in enumerate(genes) if item in ["", "."] ]
            
            # for the missing symbol positions, if the current list of symbols
            # contains a value, swap that value into the gene list.
            for x in blanks:
                try:
                    if symbols[x] != "":
                        genes[x] = symbols[x]
                except IndexError as error:
                    # very occasionally we get a variant that raises this error.
                    # The one example I've seen had ENSR=,.|ENSR00000215586.
                    # This is first split by ',' (giving ['', '.|ENSR00000215586']),
                    # then an entry is selected and split by '|'. The issue is
                    # that only the seond entry contains '|', so the lengths are
                    # discrepant
                    continue
        
        if genes is not None:
            genes = [ x if x not in ["", "."] else None for x in genes ]
        
        return genes
    
    def get_genes(self):
        """ split a gene string into list of gene names
        
        Returns:
            list of gene ID lists, one per allele
        """
        
        if self.genes is None:
            return []
        
        return self.genes
    
    def get_consequences(self, info, alts, masked):
        """ get a list of consequences for the different alt alleles
        
        Args:
            info: dictionary of keys and values for the info fields. Contains
                a 'CQ' entry, which is a comma-separated (for different alt
                alleles) and pipe-separated (for different genes) VEP-annotated
                consequence.
            alts: list of alternative alleles for the variant.
            masked: list of alternative alleles that we don't consider. These
                are identified as alt alleles with zero depth in the individual.
                This can occur due to multi-sample calling.
        
        Returns:
            list of gene lists, one per alternative allele (after removing the
            masked alt alleles.)
        """
        
        pos = [ i for i, x in enumerate(alts) if x not in masked ]
        cq = None
        if "CQ" in info:
            cq = info["CQ"].split(',')
            cq = [ cq[i].split('|') for i in pos ]
        
        # Allow for sites at the end of exons, changing from a conserved base.
        # These haven't been annotated in the VCF, so we modify the VEP
        # consequence. We only need to account for "missense_variant" and
        # "splice_region_variant", since variants near an exon end can only have
        # these consequences. The splice_region variants would ordinarily be
        # missed. We might erroneously change missense_variants in transcripts
        # where in one transcript the exon ends, while the other transcript the
        # exon continues, but those seem sufficiently rare.
        if (self.get_chrom(), self.get_position()) in self.last_base:
            required = ["missense_variant", "splice_region_variant"]
            new = "conserved_exon_terminus_variant"
            
            for old in required:
                cq = [ [ z.replace(old, new) for z in x ] for x in cq ]
        
        return cq
    
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
            return [ l for sublist in self.consequence for l in sublist ]
        
        # At one point, the VCFs lacked per gene consequences, but could have
        # multiple gene symbols (if they lacked a HGNC field but did have a
        # HGNC_ALL field). These variants will have multiple genes, but only one
        # consequence. Return the consequence as is, in order to retain the
        # same output for those VCFs.
        if len(self.get_genes()) > 1 and len(self.consequence) == 1:
            return self.consequence[0]
        
        # find the consequernce terms for the given HGNC symbol. The HGNC
        # symbols and consequences are lists of symbols/consequences per allele.
        # We need to look in the gene lists to first to identify the allele
        # index position, then the nested symbol position, so we can extract the
        # correct consequence term.
        cq = []
        for x, item in enumerate(self.get_genes()):
            if hgnc_symbol not in item:
                continue
            
            cq.append(self.consequence[x][item.index(hgnc_symbol)])
        
        return cq
    
    def get_zero_count_alleles(self, info, alt_alleles):
        ''' get a list of alleles with zero counts
        
        Some variants have multiple alts, so we need to select the alt with
        the most severe consequence. However, in at least one version of the
        VCFs, one of the alts could have zero counts, which I believe resulted
        from the population based multi-sample calling. We need to drop the
        consequences recorded for zero-count alternate alleles before finding
        the most severe.
        
        Args:
            info: info dictionary, which may contain an 'AC' key, where the
                values are a comma-separated list of counts for the alternate
                alleles. Some variants lack this field (such as CNVs).
            alt_alleles: tuple of alt alleles
        
        Returns:
            list of alleles with zero depth
        '''
        
        if 'AC' in info:
            counts = info['AC'].split(',')
            assert len(counts) == len(alt_alleles)
            
            # find the positions of alleles where the allele count is zero
            pos = [ i for i, x in enumerate(counts) if x == '0' ]
            
            # return the alleles with zero-count ,so we can mask them out
            return [ alt_alleles[i] for i in pos ]
        
        return []
    
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
        
        if self.mnv_code is not None:
            if self.mnv_code == 'masked_stop_gain_mnv':
                cq = [ x for x in cq if x != 'stop_gained' ]
            elif self.mnv_code == 'modified_stop_gained_mnv':
                cq.append('stop_gained')
        
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
        
        if self.mnv_code is not None:
            if self.mnv_code == 'modified_synonymous_mnv':
                cq = []
            elif self.mnv_code == 'modified_protein_altering_mnv':
                cq.append('missense_variant')
            elif self.mnv_code == 'masked_stop_gain_mnv':
                cq.append('missense_variant')
        
        return len(set(cq) & self.missense_consequences) > 0
    
    def is_synonymous(self, hgnc_symbol=None):
        """ checks if a variant has a missense-styled consequence
        """
        
        if self.consequence is None:
            return False
        
        cq = self.get_per_gene_consequence(hgnc_symbol)
        
        return len(set(cq) - self.synonymous_consequences) == 0
    
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
    
    
