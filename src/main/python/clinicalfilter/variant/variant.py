""" classes for handling variants
"""

class Variant(object):
    """ generic functions for variants
    """
    
    # define some codes used in ped files to identify male and female sexes
    male_codes = set(["1", "m", "M", "male"])
    female_codes = set(["2", "f", "F", "female"])
    
    x_pseudoautosomal_regions = [(60001, 2699520), (154930290, 155260560), \
        (88456802, 92375509)]
    y_pseudoautosomal_regions = [(10001, 2649520), (59034050, 59363566)]
    
    def __init__(self, chrom, position, variant_id, ref_allele, alt_allele, filt):
        """ initialise the object with the definition values
        """
        
        self.chrom = chrom
        self.position = int(position)
        
        self.variant_id = variant_id
        self.mutation_id = "NA"
        self.set_mutation_id(variant_id)
        
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        
        # self.quality = quality
        self.filter = filt
        
        # intialise variables that will be set later
        self.gender = None
        self.vcf_line = None
        self.format = None
        self.inheritance_type = None
        self.info = {}
        self.genes = None
        self.consequence = None
        self.genotype = None
        
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
    
    def __str__(self):
        
        return "{0} chr{1}: {2} {3} in {4}".format(str(self.__class__.__name__), \
                  self.chrom, self.position, self.genotype, self.genes)
    
    def add_format(self, format_keys, sample_values):
        """Parses the FORMAT column from VCF files.
        
        Args:
            format_keys: FORMAT text from a line in a VCF file
            sample_values: the values for the format keys
        """
        
        self.format = {}
        
        tag_labels = format_keys.split(":")
        tag_values = sample_values.split(":")
        
        for i, value in enumerate(tag_values):
            self.format[tag_labels[i]] = value
    
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
