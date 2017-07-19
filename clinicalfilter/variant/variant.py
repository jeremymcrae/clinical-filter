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

from clinicalfilter.variant.info import Info

class Variant(object):
    """ generic functions for variants
    """
    
    # define some codes used in ped files to identify male and female sexes
    male_codes = set(["1", "m", "M", "male"])
    female_codes = set(["2", "f", "F", "female"])
    
    x_pseudoautosomal_regions = [(60001, 2699520), (154930290, 155260560), \
        (88456802, 92375509)]
    y_pseudoautosomal_regions = [(10001, 2649520), (59034050, 59363566)]
    known_genes = None
    
    @classmethod
    def set_known_genes(cls_obj, known_genes):
        cls_obj.known_genes = known_genes
    
    def __init__(self, chrom, position, id, ref, alts, filter, info=None,
            format=None, sample=None, gender=None, mnv_code=None):
        """ initialise the object with the definition values
        """
        
        self.chrom = chrom
        self.position = int(position)
        
        self.variant_id = id
        self.mutation_id = "NA"
        self.set_mutation_id(self.variant_id)
        
        self.ref_allele = ref
        self.alt_alleles = tuple(alts.split(','))
        
        self.mnv_code = mnv_code
        self.filter = filter
        
        # intialise variables that will be set later
        self.genes = None
        self.consequence = None
        self.inheritance_type = None
        
        self.gender = None
        if gender is not None:
            self._set_gender(gender)
        
        self.vcf_line = None
        
        self.format = None
        if format is not None and sample is not None:
            self.add_format(format, sample)
        
        self.info = Info(info, self.mnv_code)
        masked = self.get_low_depth_alleles(self.ref_allele, self.alt_alleles)
        self.set_genes_and_consequence(masked, self.get_chrom(),
            self.get_position(), self.ref_allele, self.alt_alleles)
        
        self.genotype = None
        if self.format is not None and self._get_gender() is not None:
            self.set_genotype()
    
    def is_lof(self, gene_symbol=None):
        return self.info.is_lof(gene_symbol)
    def is_missense(self, is_cnv, gene_symbol=None):
        return self.info.is_missense(is_cnv, gene_symbol)
    def is_synoymous(self, gene_symbol=None):
        return self.info.is_synoymous(gene_symbol)
    
    def __repr__(self):
        ''' repr function for Variant objects. SNV(...) and CNV(...) also work
        '''
        
        def quote(value):
            if value is not None:
                value = '"{}"'.format(value)
            return value
        
        # reprocess the format dictionary back to the original text strings
        keys, sample = None, None
        if self.format is not None:
            keys = quote(':'.join(sorted(self.format)))
            sample = quote(':'.join([ self.format[x] for x in sorted(self.format) ]))
        
        info = quote(self.info)
        gender = quote(self.gender)
        mnv_code = quote(self.mnv_code)
        
        return '{}(chrom="{}", position={}, id="{}", ref="{}", alts="{}", ' \
            'filter="{}", info={}, format={}, sample={}, gender={}, ' \
            'mnv_code={})'.format(type(self).__name__, self.chrom,
            self.position, self.variant_id, self.ref_allele,
            ','.join(self.alt_alleles), self.filter, info, keys, sample,
            gender, mnv_code)
    
    def __hash__(self):
        return hash(str(self))
    
    def __eq__(self, other):
        return hash(self) == hash(other)
    
    def _set_gender(self, gender):
        """ sets the gender of the individual for the variant
        """
        
        if gender in self.male_codes:
            self.gender = "male"
        elif gender in self.female_codes:
            self.gender = "female"
        else:
            raise ValueError("unknown gender code")
        
        self.set_inheritance_type(self.get_position(), self.is_male())
    
    def _get_gender(self):
        """returns the gender for a person (1, M = male, 2, F = female).
        """
        return self.gender
    
    def set_mutation_id(self, variant_id):
        """ sets the mutation ID based on the VCF ID field
        
        The variant ID can be either "." for null value, an rsID, a HGMD ID,
        a COSMIC ID, or any combination of those (including multiple HGMD IDs
        for a single variant).
        
        Args:
            variant_id: string from the VCF ID field, can be rsID, or a list of
                &-separated IDs, which can include COSMIC and HGMD IDs.
        """
        
        if variant_id != ".":
            variant_id = variant_id.split("&")
            ids = []
            for value in variant_id:
                # include everything that isn't an rsID
                if not value.startswith("rs"):
                    ids.append(value)
            
            if len(ids) > 0:
                self.mutation_id = ",".join(ids)
                    
    def get_mutation_id(self):
        return self.mutation_id
    
    def is_male(self):
        """ returns True/False for whether the person is male
        """
        
        return self._get_gender() in self.male_codes
    
    def add_format(self, keys, values):
        """Parses the FORMAT column from VCF files.
        
        Args:
            keys: FORMAT text from a line in a VCF file
            values: the values for the format keys
        """
        
        self.format = dict(zip(keys.split(":"), values.split(":")))
    
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
        
        allele_counts = ['1'] * len(alts)
        if 'AC' in self.info:
            allele_counts = self.info['AC'].split(',')
        
        allele_depths = ['10'] * len(alts)
        if 'AD' in self.format:
            allele_depths = self.format['AD'].split(',')[1:]
        
        counts = list(zip(allele_counts, allele_depths))
        
        assert len(counts) == len(alts)
        
        # find the positions of alleles where the allele count is zero,
        # or indels with 1 alt read
        pos = set()
        for i, (count, depth) in enumerate(counts):
            if count == '0':
                pos.add(i)
            elif depth == '1' and is_indel(ref, alts[i]):
                pos.add(i)
        
        # return the alleles with zero-count ,so we can mask them out
        return [ alts[i] for i in sorted(pos) ]
    
    def add_vcf_line(self, vcf_line):
        self.vcf_line = vcf_line
    
    def get_vcf_line(self):
        return self.vcf_line
        
    def set_inheritance_type(self, pos, is_male):
        """ sets the chromosome type (eg autosomal, or X chromosome type).
        
        provides the chromosome type for a chromosome (eg Autosomal, or
        X-chrom male etc). This only does simple string matching. The
        chromosome string is either the chromosome number, or in the case of
        the sex-chromosomes, the chromosome character. This doesn't allow for
        chromosomes to be specified as "chr1", and sex chromosomes have to be
        specified as "X" or "Y", not "23" or "24".
        
        Args:
            pos: position on the chromosome
            is_male: True/False for whether the individual is male
        """
        
        if self.get_chrom() not in ["chrX", "ChrX", "X", "chrY", "ChrY", "Y"]:
            self.inheritance_type = "autosomal"
        elif self.get_chrom() in ["chrX", "ChrX", "X"]:
            # check if the gene lies within a pseudoautosomal region
            for start, end in self.x_pseudoautosomal_regions:
                if start < pos < end:
                    self.inheritance_type = "autosomal"
                    return
            
            if is_male:
                self.inheritance_type =  "XChrMale"
            else:
                self.inheritance_type = "XChrFemale"
        elif self.get_chrom() in ["chrY", "ChrY", "Y"]:
            # check if the gene lies within a pseudoautosomal region
            for start, end in self.y_pseudoautosomal_regions:
                if start < pos < end:
                    self.inheritance_type = "autosomal"
                    return
            if is_male:
                self.inheritance_type =  "YChrMale"
            else:
                self.inheritance_type = "YChrFemale"
    
    def get_inheritance_type(self):
        """ return the variant chromosomal inheritance type
        """
        
        return self.inheritance_type
    
    def get_chrom(self):
        """ return the variant chromosome
        """
        
        return self.chrom
    
    def get_position(self):
        """ return the variant chromosomal position
        """
        
        return self.position
    
    def get_genotype(self):
        """ return the genotype value
        """
        
        return self.genotype
