""" loads variant data from individuals into a single object, so we can easily
get genotype data for a single variant from all the family members.
"""

from clinicalfilter.variant_cnv import CNV

class TrioGenotypes(object):
    """ a class to hold genotypes for the members of a trio
    """

    def __init__(self, child_variant):
        """ initiate the class with the childs variant

        Args:
            child_variant: Variant object
        """

        self.child = child_variant

        self.chrom = self.child.get_chrom()
        self.position = self.child.get_position()
        self.inheritance_type = self.child.inheritance_type
        self.gene = self.child.gene

    def convert_chrom_to_int(self, chrom):
        """ converts a chromosome string to an int (if possible) for sorting.

        Args:
            chrom: string (eg "1", "2", "3" ... "22", "X", "Y")

        Returns:
            int value of chrom
        """

        # set the integer values of the sex chromosomes
        chrom_dict = {"X": 23, "CHRX": 23, "Y": 24, "CHRY": 24, "MT": 25, \
            "CHRMT": 25}

        try:
            chrom = int(chrom)
        except ValueError:
            chrom = chrom_dict[chrom.upper()]

        return chrom

    def __eq__(self, other):
        return (self.convert_chrom_to_int(self.chrom), int(self.position)) == \
              (self.convert_chrom_to_int(other.chrom), int(other.position))

    def __ne__(self, other):
        return (self.convert_chrom_to_int(self.chrom), int(self.position)) != \
              (self.convert_chrom_to_int(other.chrom), int(other.position))

    def __lt__(self, other):
        return (self.convert_chrom_to_int(self.chrom), int(self.position)) < \
              (self.convert_chrom_to_int(other.chrom), int(other.position))

    def __repr__(self):
        return self.__str__()

    def __str__(self):

        chrom = self.get_chrom()
        position = self.get_position()
        (child, mother, father) = self.get_trio_genotype()

        return "chr%s: %s - %s%s%s" % (chrom, position, child, mother, father)

    def __hash__(self):
        return hash(self.__repr__())

    def is_cnv(self):
        """ checks whether the variant is for a CNV
        """

        return isinstance(self.child, CNV)

    def add_father_variant(self, father_variant):
        self.father = father_variant

    def add_mother_variant(self, mother_variant):
        self.mother = mother_variant

    def get_position(self):
        return self.position

    def get_inheritance_type(self):
        return self.inheritance_type

    def get_trio_genotype(self):
        child_geno = self.child.get_genotype()

        if hasattr(self, "mother"):
            mother_geno = self.mother.get_genotype()
        else:
            mother_geno = "NA"

        if hasattr(self, "father"):
            father_geno = self.father.get_genotype()
        else:
            father_geno = "NA"

        return (child_geno, mother_geno, father_geno)

    def get_chrom(self):
        return self.chrom

    def get_gene(self):
        return self.gene

    def get_de_novo_genotype(self):
        """ get the de novo genotype combination for the chromosome/sex
        """

        # set the standard de novo genotype combination
        de_novo_genotype = (1, 0, 0)
        
        # account for X chrom de novos in males
        if self.inheritance_type == "XChrMale":
            de_novo_genotype = (2, 0, 0)
        
        return de_novo_genotype
    
    def passes_de_novo_checks(self, pp_filter):
        """ checks if the child's de novo variants passes filters

        Some variants are de novo in the child, and if they are, then we should
        subject them to additional filtering to see if they have been passed by
        denovogear. Note that both parents have genotypes specified, as we
        have checked at an earlier stage that the parents are present.
        
        Args:
            pp_filter: float between 0 and 1, being the threshold for the PP_DNM filter

        Returns:
            boolean value for whether the variant should be included
        """

        # currently hard code the filtering fields. The de novo field indicates
        # whether the variant is de novo, I don't know how this is assigned. The
        # project filter field indicates an internal filter, curently whether
        # the variant passed MAF, alternate frequency, and segmental duplication
        # criteria.
        de_novo_field = set(["DENOVO-SNP", "DENOVO-INDEL"])

        # if the variant is not de novo, don't worry about de novo filtering
        if self.get_trio_genotype() != self.get_de_novo_genotype():
            return True

        # check the VCF record to see whether the variant has been screened out.
        # Either DENOVO-SNP or DENOVO-INDEL should be in the info.
        if len(set(self.child.info) & de_novo_field) < 1:
            return False

        if "PP_DNM" in self.child.format and \
                float(self.child.format["PP_DNM"]) < pp_filter:
            return False

        return True



