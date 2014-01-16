""" classes for handling variants
"""

class Variant(object):
    """ generic functions for variants
    """
    
    def __init__(self, chrom, position, snp_id, ref_allele, alt_allele, quality, filter):
        """ initialise the object with the definition values
        """
        
        self.chrom = chrom
        self.position = position
        self.id = snp_id
        self.set_mutation_id()
        
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        
        self.quality = quality
        self.filter = filter
        
        # define some codes used in ped files to identify male and female sexes
        self.male_codes = set(["1", "m", "M", "male"])
        self.female_codes = set(["2", "f", "F", "female"])
        
        self.pseudoautosomal_regions = [(1,2699520), (154930290,155260560), (88456802,92375509)]
        
    def set_gender(self, gender):
        """ sets the gender of the individual for the variant
        """
        if gender in self.male_codes:
            self.gender = "male"
        elif gender in self.female_codes:
            self.gender = "female" 
        else:
            raise ValueError("unknown gender")
        
        self.set_inheritance_type()
    
    def set_mutation_id(self):
        """ sets the mutation ID based on the ID column value
        """
        
        # the variant ID can be either "." for null value, an rsID, a HGMD ID,
        # a COSMIC ID, or any combination of those (including multiple HGMD IDs
        # for a single variant). 
        mutation_id = self.id
        
        if mutation_id == ".":
            self.mutation_id = "NA"
        else:
            mutation_id = mutation_id.split("&")
            ids = []
            for value in mutation_id:
                # include everything that isn't an rsID
                if not value.startswith("rs"):
                    ids.append(value)
            
            if len(ids) == 0:
                self.mutation_id = "NA"
            else:
                self.mutation_id = ",".join(ids)
                    
    def get_mutation_id(self):
        return self.mutation_id
    
    def get_gender(self):
        """returns the gender for a person (1, M = male, 2, F = female).
        """
        return self.gender
    
    def is_male(self):
        """ returns True/False for whether the person is male
        """
        
        return self.get_gender() in self.male_codes
    
    def is_female(self):
        """ returns True/False for whether the person is male
        """
        
        return self.get_gender() in self.female_codes
    
    def __str__(self):
        
        string = "%s chr%s: %s %s in %s" % (str(self.__class__.__name__), self.chrom, \
                 self.position, self.genotype, self.gene)
        return string
    
    def add_format(self, format_keys, sample_values):
        """Parses the FORMAT column from VCF files.
        
        Args:
            format_keys: FORMAT text from a line in a VCF file
            sample_values: the values for the format keys
        """
        
        self.format = {}
        
        tag_labels = format_keys.split(":") # the first item in formats are the tags DP:FP:ETC
        tag_values = sample_values.split(":") # the second item are the corresponding values.
        
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
        
        if self.chrom not in ["chrX", "ChrX", "X"]:
            self.inheritance_type = "autosomal"
        else:
            # check if the gene lies within a pseudoautosomal region
            for start, end in self.pseudoautosomal_regions:
                if start < int(self.position) < end or start < int(self.position) < end:
                    self.inheritance_type = "autosomal"
            
            if self.is_male():
                self.inheritance_type =  "XChrMale"
            if self.is_female():
                self.inheritance_type = "XChrFemale"
    
    def get_inheritance_type():
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
        
        return str(self.position)
    
    def get_genotype(self):
        """ return the genotype value
        """
        
        return self.genotype
    
    def find_max_allele_frequency(self, populations):
        """gets the maximum allele frequency for a variant in a VCF record
        
        Finds the maximum allele frequency recorded for a variant across
        different populations.
        
        Args:
            populations: list of population IDs to search
          
        Returns:
            the maximum allele frequency found within the populations in the
            variant record
        """
        
        max_allele_frequency = -100
        # run through all the possible populations in the VCF file (typically the 1000 Genomes 
        # populations (AFR_AF, EUR_AF etc), an internal popuation (DDD_AF), and a AF_MAX field)
        for key in populations:
            if key in self.info:
                number = self.get_number(self.info[key])
                if not self.is_number(number):
                    continue
                # if number > 0.5:
                #     number = 1 - number
                #     record[key] = str(number)
                if number > max_allele_frequency:
                    max_allele_frequency = number
        
        # return NA for variants without MAF recorded
        if max_allele_frequency == -100:
            max_allele_frequency = "NA"
        
        return str(max_allele_frequency)

