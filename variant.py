""" simple class for handling biallelic variants
"""

import math

class VcfInfo(object):
    """ parses the VCF info field, and checks whether the variant passes 
    filtering criteria.
    """
    
    def __init__(self):
        """
        """
        
        self.show_fail_point = False
    
    def add_info(self, info_values):
        """Parses the INFO column from VCF files.
        
        Args:
            info_values: INFO text from a line in a VCF file
        """
        
        self.info = {}
        
        for item in info_values.split(";"):
            if "=" in item:
                key, value = item.split("=")
            else:
                key, value = item, True
            self.info[key] = value
        
        # add the filter value, as we filter with the info dict
        self.info["FILTER"] = self.filter
        
        # make sure we have gene and consequence keys in the info dict, for 
        # the filter to work with
        if "HGNC" not in self.info:
            self.info["HGNC"] = None
        
        self.gene = self.info["HGNC"]
        
        if "CQ" not in self.info:
            self.info["CQ"] = None
    
    def get_number(self, values):
        """ converts a string into a number
        
        This function is used for GAPI VCF files which might have multiple files
        seprated by "," at INFO columns. This should be used for generic VCF
        as the end user may not know the only the first value is returned 
        """
        # if the string can be directly converted to a float, simply return that
        try:
            value = float(values)
        # occasionally we get comma-separated pairs (eg '.,0.639860'). Try to convert each of these 
        # in turn, if any can be converted to floats, return that value
        except ValueError:
            values = values.split(",")
            for value in values:
                try:
                    value = float(value)
                    break
                except ValueError:
                    pass
        except:
            value = values
        return value
    
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
    
    def show_fail(self, key, value, condition, filter_values):
        print(str(key) + ": " + str(value) + " not " + str(condition) + \
                  " " + str(filter_values))
    
    def passes_filters(self, filters):
        """Checks whether a VCF record passes user defined criteria.
        
        Args:
            record: A dictionary entry for a single variant converted from a
            VCF file.
            
        Returns:
            boolean value for whether the record passes the filters
        """
        
        self.show_fail_point = False
        if self.chrom == '1' and self.position == '18937':
            show_fail_point = True
        
        passes = True
        for key in self.info:
            if key not in filters:
                continue
            
            value = self.info[key]
            condition = filters[key][0]
            filter_values = filters[key][1]
            
            if condition == "list" and self.fails_list(value, filter_values):
                passes = False
                break
            elif condition == "greater_than" and self.fails_greater_than(value, filter_values):
                passes = False
                break
            elif condition == "smaller_than" and self.fails_smaller_than(value, filter_values):
                passes = False
                break
            elif condition == "equal" and self.fails_equal(value, filter_values):
                passes = False
                break
            elif condition == "not" and self.fails_not(value, filter_values):
                passes = False
                break
            elif condition == "multiple_not" and self.fails_multiple_not(value, filter_values):
                passes = False
                break
            elif condition == "startswith" and self.fails_start_string(value, filter_values):
                passes = False
                break
            elif condition == "endswith" and self.fails_end_string(value, filter_values):
                passes = False
                break
            elif condition == "range" and self.fails_range(value, filter_values):
                passes = False
                break
                
            if passes == False:
                break
        
        if passes == False:
            if self.show_fail_point:
                self.show_fail(key, value, condition, filter_values)
        
        return passes
    
    def fails_list(self, value, filter_values):
        """ checks whether the vcf value is within a list 
        """
        
        fails = False
        if value not in filter_values:
            fails = True
        
        return fails
    
    def fails_greater_than(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        value = self.get_number(value)
        try:
            value > filter_values
        except TypeError:
            return False
            
        fails = False
        if value < filter_values and self.is_number(value):
            fails = True
        
        return fails
    
    def fails_smaller_than(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        # some of the MAF values are 1 - MAF due to being for a population that was 
        # genotyped on the opposing strand. We need to convert those back.
        # if key in self.tags_dict["MAX_MAF"]:
        #     value = self.get_number(value)
        #     if self.is_number(value):
        #         if value > 0.5:
        #             value = 1 - value
        value = self.get_number(value)
        try:
            value > filter_values
        except TypeError:
            return False
        
        fails = False
        if value > filter_values and self.is_number(value):
            fails = True
        
        return fails
                    
    def fails_equal(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        fails = False
        if value != filter_values:
            fails = True
            
        return fails
    
    def fails_not(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        fails = False
        if value == filter_values:
            fails = True
            
        return fails
    
    def fails_range(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        fails = False
        start, end = filter_values
        if self.get_number(value) < start or self.get_number(value) > end:
            fails = True
            
        return fails
    
    def fails_start_string(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        fails = False
        if not value.startswith(filter_values):
            fails = True
        
        return fails
        
    def fails_end_string(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        fails = False
        if not value.endswith(filter_values):
            fails = True
        
        return fails
        
    def fails_multiple_not(self, value, filter_values):
        """ checks whether values are not within a filter range
        """
        
        # sometimes we want to make sure two related keys don't contain certain
        # values. Specifically, we want to exclude variants where polyphen 
        # annotation is "benign" and sift annotation is "tolerated". These are
        # provided as a list of  tuples. We catch the first key, and then check
        # the second key at the same time
        has_all_values = True
        for key, filter_value in filter_values:
            # pull out the value for each key in the list, and check whether it
            # contains the filter value
            temp_value = self.info[key]
            if temp_value is None:
                has_all_values = False
            # the value should be something like benign(0.05), or tolerated(0.3),
            # and we just check if "benign", or "tolerated" are in the 
            # corresponding values
            elif filter_value not in temp_value:
                has_all_values = False
        
        fails = False
        if has_all_values:
            passes = True
        
        return fails


class Variant(VcfInfo):
    """ generic functions for variants
    """
    
    def __init__(self, chrom, position, snp_id, ref_allele, alt_allele, quality, filter):
        """ initialise the object with the definition values
        """
        
        self.chrom = chrom
        self.position = position
        self.id = snp_id
        
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
        chromosomes to be specified as 'chr1', and sex chromosomes have to be 
        specified as 'X' or 'Y', not '23' or '24'.
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
    
    def get_position():
        """ return the variant chromosomal inheritance type
        """
        
        return self.position
    
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
    
class SNV(Variant):
    """ a class to take a genotype coded as 0, 1 or 2, and be able to perform simple functions, like
    reporting whether it is heterozygous, homozygous, or neither, depending on whether the variant 
    is on the X chromosome, and if so, whether the individual is male or female.
    """
    
    def get_key(self):
        """ return a tuple to identify the variant
        """
        
        return (self.chrom, self.position)
    
    def set_genotype(self, genotype=None):
        """ sets the genotype of the variant
        """
        
        if hasattr(self, "format"):
            self.genotype = self.format["GT"]
            self.genotype = self.convert_genotype()
        elif genotype is not None:
            self.genotype = genotype
        else:
            raise NotImplementedError("cannot find a genotype")
        
        self.set_reference_genotypes()
        self.convert_genotype_code_to_alleles()
    
    def set_default_genotype(self):
        """ for variants lacking genotypes, set a defaul genotype
        """
        
        self.genotype = 0
        
        self.set_reference_genotypes()
        self.convert_genotype_code_to_alleles()
    
    def convert_genotype(self):
        """Maps genotypes from two character format to single character.
        
        Args:
            GT: genotype in two character format. eg '0/0'
        
        Returns:
            Genotype in single character format. eg '0'
        """
        # This function might run quicker as the following. Possibly could speed up further by passing 
        # reference to a dictionary, rather than creating new each time.
        genotype_dict = {'00': 0, '01': 1, '10': 1, '12': 1, '21': 1, '02': 1, '20': 1, '11': 2, '22': 2}
        return genotype_dict[self.genotype[0] + self.genotype[-1]]
    
    def convert_genotype_code_to_alleles(self):
        """ converts a 0/1/2 genotype code to a set of alleles, depending on the chromosome
        """
        
        if self.inheritance_type == "autosomal":
            self.convert_autosomal_genotype_code_to_alleles()
        elif self.inheritance_type == "XChrMale" or self.inheritance_type == "XChrFemale":
            self.convert_allosomal_genotype_code_to_alleles()
    
    def set_reference_genotypes(self):
        """ sets reference genotypes for homozygotes and heterozygotes, for checking against
        """
        
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale':
            self.hom_ref = set([self.ref_allele, self.ref_allele])
            self.het = set([self.ref_allele, self.alt_allele])
            self.hom_alt = set([self.alt_allele, self.alt_allele])
        elif self.inheritance_type == 'XChrMale':
            self.hom_alt = set([self.alt_allele])
            self.hom_ref = set([self.ref_allele])
    
    def is_het(self):
        """ returns whether a variant is heterozygous
        """
        
        is_het = False
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale':
            if self.alleles == self.het:
                is_het = True
        elif self.inheritance_type == 'XChrMale':
            pass
        
        return is_het
    
    def is_hom_alt(self):
        """ returns whether a genotype is homozygous for the alternate allele
        """
        
        is_hom_alt = False
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale' or self.inheritance_type == 'XChrMale':
            if self.alleles == self.hom_alt:
                is_hom_alt = True
        
        return is_hom_alt
    
    def is_hom_ref(self):
        """ returns whether a variant is homozygous for the reference allele
        """
        
        is_hom_ref = False
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale' or self.inheritance_type == 'XChrMale':
            if self.alleles == self.hom_ref:
                is_hom_ref = True
        
        return is_hom_ref
    
    def is_not_ref(self):
        """ returns whether a variant is not the homozygous for the reference allele
        """
        
        is_not_hom_ref = False
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale' or self.inheritance_type == 'XChrMale':
            if self.alleles != self.hom_ref:
                is_not_hom_ref = True
        
        return is_not_hom_ref
    
    def is_not_alt(self):
        """ returns whether a variant is not homozygous for the alternate allele
        """
        
        is_not_hom_alt = False
        if self.inheritance_type == 'autosomal' or self.inheritance_type == 'XChrFemale' or self.inheritance_type == 'XChrMale':
            if self.alleles != self.hom_alt:
                is_not_hom_alt = True
        
        return is_not_hom_alt
    
    def convert_autosomal_genotype_code_to_alleles(self):
        """converts a genotype code to a set of alleles
        
        returns:
            an unsorted set of allele codes
        """
        
        genotype = str(self.genotype)
        
        if genotype == '0':
            self.alleles = set([self.ref_allele, self.ref_allele])
        elif genotype == '1':
            self.alleles = set([self.ref_allele, self.alt_allele])
        elif genotype == "2":
            self.alleles = set([self.alt_allele, self.alt_allele])
        else:
            raise ValueError("genotype code '" + genotype + "' is not recognised")
        
    def convert_allosomal_genotype_code_to_alleles(self):
        """converts a genotype on the x-chromosome into the set of alleles
        
        returns:
            an unsorted set of allele codes
        """
        
        genotype = str(self.genotype)
        
        if self.is_male():
            if genotype == "0":
                self.alleles = set([self.ref_allele])
            elif genotype == "2":
                self.alleles = set([self.alt_allele])
            elif genotype == "1":
                self.alleles = set([self.ref_allele, self.alt_allele])
                raise ValueError("heterozygous X-chromomosome male")
            else:
                raise ValueError("unknown genotype '" + str(genotype))
        elif self.is_female():
            self.convert_autosomal_genotype_code_to_alleles()
        else:
            raise ValueError("Unknown gender: " + self.gender)

class CNV(Variant):
    
    def set_genotype(self):
        """ sets the genotype of the variant
        """
        
        if self.alt_allele == "<DUP>":
            self.genotype = "dup"
        elif self.alt_allele == "<DEL>":
            self.genotype = "del"
        elif self.alt_allele == "no_variation":
            self.alt_allele = "ref"
        
        self.set_range()
        self.set_reference_genotypes()
    
    def set_default_genotype(self):
        """ set a default genotype for individuals without one
        """
        
        self.genotype = "ref"
        self.set_range()
        self.set_reference_genotypes()
    
    def set_reference_genotypes(self):
        """ sets reference genotypes for homozygotes and heterozygotes, for checking against
        """
        
        self.ref_genotypes = set(["ref"])
        self.alt_genotypes = set(["del", "dup"])
    
    def set_range(self):
        """ sets the range for the CNV
        """
        
        self.start_position = self.position
        
        if hasattr(self, "info"):
            self.end_position = self.info["END"]
        else:
            self.end_position = str(int(self.start_position) + 10000)
        
        self.range = (self.start_position, self.end_position)
        self.size = int(self.end_position) - int(self.start_position)
    
    def get_key(self):
        """ return a tuple to identify the variant
        """
        
        try:
            return (self.chrom, self.start_position, self.end_position)
        except AttributeError:
            self.set_range()
        
        return (self.chrom, self.start_position, self.end_position)
    
    def get_genotype(self):
        """ return the genotype value
        """
        
        return self.genotype
    
    def calculate_cnv_size_tolerance(self):
        """ calculates the size range of CNVs that might match a given CNV size.
        
        Args:
            cnv_size: int value for CNV size in base pairs
        
        Returns:
            tuple of minimum and maximum sizes in base pairs
        """
        
        min_size = self.size - abs(100 * math.sqrt(self.size + 2500)) + 5000
        max_size = self.size + 100 * math.sqrt(self.size)
        
        return (min_size, max_size)
    
    def is_het(self):
        return False
        # if self.genotype in self.alt_genotypes:
        #     return True
        # else:
        #     return False
    
    def is_hom_alt(self):
        return False
        # if self.genotype in self.alt_genotypes:
        #     return True
        # else:
        #     return False
    
    def is_hom_ref(self):
        return True
        # if self.genotype in self.ref_genotypes:
        #     return True
        # else:
        #     return False
    
    def is_not_ref(self):
        return False
        # if self.genotype in self.alt_genotypes:
        #     return True
        # else:
        #     return False
    
    def is_not_alt(self):
        return True
        # if self.genotype in self.ref_genotypes:
        #     return True
        # else:
        #     return False

