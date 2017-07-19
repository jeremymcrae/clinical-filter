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

from clinicalfilter.variant.symbols import Symbols

class Info(object):
    """ parses the VCF INFO field
    """
    
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
    last_base = set([])
    populations = []
    
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
        
        self.mnv_code = mnv_code
        self.info = {}
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
    
    def set_genes_and_consequence(self, masked):
        ''' find the gene symbols and consequences for good alleles
        '''
        self.genes = self.get_gene_from_info(self.info, self.alt_alleles, masked)
        self.consequence = self.get_consequences(self.info, self.alt_alleles, masked)
    
    def __str__(self):
        ''' reprocess the info dictionary back into a string, correctly sorted
        '''
        
        info = []
        for key, value in sorted(self.info.items()):
            entry = key
            if value != True:
                entry = '{}={}'.format(key, value)
            info.append(entry)
            
        return ';'.join(info)
    
    def __getitem__(self, key):
        return self.info[key]
    
    def __setitem__(self, key, value):
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
    
    def __contains__(self, key):
        return key in self.info
    
    def parse_gene_symbols(self, alts, masked):
        """ parses the available gene symbols in the INFO.
        
        Args:
            alts: list of alternative alleles for the variant
            masked: list of alternative alleles that we don't consider. These
                are identified as alt alleles with zero depth in the individual.
                This can occur due to multi-sample calling.
        
        Returns:
            list of gene lists, one per alternative allele (after removing the
            masked alt alleles.)
        """
        
        pos = [ i for i, x in enumerate(alts) if x not in masked ]
        return [ Symbols(self.info, i) for i in pos ]
    
    def get_genes(self):
        """ split a gene string into list of gene names
        
        Returns:
            list of gene ID lists, one per allele
        """
        
        if self.symbols is None:
            return []
        
        return [ x.prioritise() for x in self.symbols ]
    
    def get_consequences(self, chrom, pos, alts, masked):
        """ get a list of consequences for the different alt alleles
        
        Args:
            chrom: chromosome for the current variant
            pos: nucleotide position of the current variant
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
        if "CQ" in self.info:
            cq = self.info["CQ"].split(',')
            cq = [ cq[i].split('|') for i in pos ]
        
        # Allow for sites at the end of exons, changing from a conserved base.
        # These haven't been annotated in the VCF, so we modify the VEP
        # consequence. We only need to account for "missense_variant" and
        # "splice_region_variant", since variants near an exon end can only have
        # these consequences. The splice_region variants would ordinarily be
        # missed. We might erroneously change missense_variants in transcripts
        # where in one transcript the exon ends, while the other transcript the
        # exon continues, but those seem sufficiently rare.
        if (chrom, pos) in self.last_base:
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
        
        # find the consequence terms for the given HGNC symbol. The HGNC
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
    
<<<<<<< HEAD
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
=======
    def get_low_depth_alleles(self, ref, alts):
        ''' get a list of alleles with zero counts, or indels with 1 read
        
        Some variants have multiple alts, so we need to select the alt with
        the most severe consequence. However, in at least one version of the
        VCFs, one of the alts could have zero counts, which I believe resulted
        from the population based multi-sample calling. We need to drop the
        consequences recorded for zero-count alternate alleles before finding
        the most severe.
        
        We also want to avoid indels with only one read, because these are
        universally bad calls.
        
        Args:
            ref: reference allele
            alts: tuple of alt alleles
        
        Returns:
            list of alleles with sufficiently low depth
        '''
        
        is_indel = lambda x, y: len(x) > 1 or len(y) > 1
        
        if 'AC' in self:
            counts = self['AC'].split(',')
            assert len(counts) == len(alts)
            
            # find the positions of alleles where the allele count is zero,
            # or indels with 1 alt read
            pos = set()
            for i, x in enumerate(counts):
                if x == '0':
                    pos.add(i)
                elif x == '1' and is_indel(ref, alts[i]):
                    pos.add(i)
            
            # return the alleles with zero-count ,so we can mask them out
            return [ alt_alleles[i] for i in sorted(pos) ]
        
        return []
    
    def is_lof(self, gene_symbol=None):
>>>>>>> start isolating INFO into an independent object
        """ checks if a variant has a loss-of-function consequence
        
        Args:
            gene_symbol: HGNC symbol for which we wish to check VEP consequence.
                By default we check all the consequences listed for the variant.
        """
        
        if self.consequence is None:
            return False
        
        cq = self.get_per_gene_consequence(gene_symbol)
        
        if self.mnv_code is not None:
            if self.mnv_code == 'masked_stop_gain_mnv':
                cq = [ x for x in cq if x != 'stop_gained' ]
            elif self.mnv_code == 'modified_stop_gained_mnv':
                cq.append('stop_gained')
        
        return len(set(cq) & self.lof_consequences) > 0
    
    def is_missense(self, is_cnv, gene_symbol=None):
        """ checks if a variant has a missense-styled consequence
        
        Args:
            gene_symbol: HGNC symbol for which we wish to check VEP consequence.
                By default we check all the consequences listed for the variant.
        """
        
        if self.consequence is None:
            return False
        
        cq = self.get_per_gene_consequence(gene_symbol)
        
        if self.mnv_code is not None:
            if self.mnv_code == 'modified_synonymous_mnv':
                cq = []
            elif self.mnv_code == 'modified_protein_altering_mnv':
                cq.append('missense_variant')
            elif self.mnv_code == 'masked_stop_gain_mnv':
                cq.append('missense_variant')
        
        # CNVs can be problematic to assign VEP consequences to. Some CNVs are
        # annotated as 'coding_sequence_variant', a term which historically is
        # used in anomalous situations. SNVs no longer have a problem with this.
        missense = set(self.missense_consequences)
        if is_cnv:
            missense.add('coding_sequence_variant')
        
        return len(set(cq) & missense) > 0
    
    def is_synonymous(self, gene_symbol=None):
        """ checks if a variant has a synonymous consequence
        """
        
        if self.consequence is None:
            return False
        
        cq = self.get_per_gene_consequence(gene_symbol)
        
        return not self.is_lof(gene_symbol) and not self.is_missense(gene_symbol) and \
            len(set(cq) & self.synonymous_consequences) > 0
    
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
