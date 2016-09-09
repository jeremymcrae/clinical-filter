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

class Variant(Info):
    """ generic functions for variants
    """
    
    # define some codes used in ped files to identify male and female sexes
    male_codes = set(["1", "m", "M", "male"])
    female_codes = set(["2", "f", "F", "female"])
    
    x_pseudoautosomal_regions = [(60001, 2699520), (154930290, 155260560), \
        (88456802, 92375509)]
    y_pseudoautosomal_regions = [(10001, 2649520), (59034050, 59363566)]
    
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
            self.set_gender(gender)
        
        self.vcf_line = None
        self.info = {}
        if info is not None:
            self.add_info(info)
        
        self.format = None
        if format is not None and sample is not None:
            self.add_format(format, sample)
        
        self.genotype = None
        if self.format is not None and self.gender is not None:
            self.set_genotype()
    
    def __repr__(self):
        ''' repr function for Variant objects. SNV(...) and CNV(...) also work
        '''
        
        # reprocess the format dictionary back to the original text strings
        keys, sample = None, None
        if self.format is not None:
            keys = ':'.join(sorted(self.format))
            sample = ':'.join([ self.format[x] for x in keys.split(':') ])
            keys = '"{}"'.format(keys)
            sample = '"{}"'.format(sample)
        
        info = None
        if self.info is not None:
            info = ';'.join([ '{}={}'.format(x, self.info[x]) for x in sorted(self.info) ])
            info = '"{}"'.format(info)
        
        gender = self.gender
        if gender is not None:
            gender = '"{}"'.format(gender)
        
        mnv_code = self.mnv_code
        if mnv_code is not None:
            mnv_code = '"{}"'.format(mnv_code)
        
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
    
    def set_gender(self, gender):
        """ sets the gender of the individual for the variant
        """
        
        if gender in self.male_codes:
            self.gender = "male"
        elif gender in self.female_codes:
            self.gender = "female"
        else:
            raise ValueError("unknown gender code")
        
        self.set_inheritance_type()
    
    def get_gender(self):
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
        
        return self.get_gender() in self.male_codes
    
    def is_female(self):
        """ returns True/False for whether the person is male
        """
        
        return self.get_gender() in self.female_codes
    
    def add_format(self, keys, values):
        """Parses the FORMAT column from VCF files.
        
        Args:
            keys: FORMAT text from a line in a VCF file
            values: the values for the format keys
        """
        
        self.format = dict(zip(keys.split(":"), values.split(":")))
    
    def add_vcf_line(self, vcf_line):
        self.vcf_line = vcf_line
    
    def get_vcf_line(self):
        return self.vcf_line
        
    def set_inheritance_type(self):
        """ sets the chromosome type (eg autosomal, or X chromosome type).
        
        provides the chromosome type for a chromosome (eg Autosomal, or
        X-chrom male etc). This only does simple string matching. The
        chromosome string is either the chromosome number, or in the case of
        the sex-chromosomes, the chromosome character. This doesn't allow for
        chromosomes to be specified as "chr1", and sex chromosomes have to be
        specified as "X" or "Y", not "23" or "24".
        """
        
        if self.chrom not in ["chrX", "ChrX", "X", "chrY", "ChrY", "Y"]:
            self.inheritance_type = "autosomal"
        elif self.chrom in ["chrX", "ChrX", "X"]:
            # check if the gene lies within a pseudoautosomal region
            for start, end in self.x_pseudoautosomal_regions:
                if start < self.position < end:
                    self.inheritance_type = "autosomal"
                    return
            
            if self.is_male():
                self.inheritance_type =  "XChrMale"
            elif self.is_female():
                self.inheritance_type = "XChrFemale"
        elif self.chrom in ["chrY", "ChrY", "Y"]:
            # check if the gene lies within a pseudoautosomal region
            for start, end in self.y_pseudoautosomal_regions:
                if start < self.position < end:
                    self.inheritance_type = "autosomal"
                    return
            if self.is_male():
                self.inheritance_type =  "YChrMale"
            elif self.is_female():
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
