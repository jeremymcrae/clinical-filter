""" class for reporting found variants via text or VCF files
"""

import logging
import platform
import datetime
import gzip
import importlib
import sys

class report(object):
    """A class to report candidate variants.
    """

    def clinicalFilterVersion(self):
        """Tries to obtain the version (git tag) from clinicalfilter.version.version()
        """
        clinical_filter_version = "XXX"
        try:
            version_module = importlib.import_module("clinicalfilter.version")
            clinical_filter_version = version_module.version()
        except ImportError as ierr:
            logging.info("Failed to import clinicalfilter.version. ImportError: '{0}' (Using 'XXX')".format (ierr))
        except:
            logging.warn("Uexpected error trying to import clinicalfilter.version. '{0}'".format(sys.exc_info()[0]))
        
        return clinical_filter_version

    def printFileContent(self, path):
        """ prints the text content of a file.
        """
        f = open(path, 'r')
        for line in f:
            if not line.startswith('#'):
                logging.info("#" + line.strip())
        f.close()
    
    def printSectionTitle(self, title):
        """prints a title to standard out.
        
        Args:
            title: a string naming an inheritance model, such as 'Autosomal Recessive HOM'.
            
        Returns:
            Nothing, just prints to the screen.
        """
        logging.info(75 * "*")
        logging.info(title.upper())
        logging.info(75 * "*")
    
    def printReportInputInfo(self, title):
        """report the parameters used to run the script, such as python version, script version, 
        filters and weights, as well as which columns were chosen for output. 
        """
        
        python_version = "%s %s %s" % (platform.python_version(), platform.python_build(), platform.python_compiler())
        
        # capture the program title
        logging.info(75 * "#")
        logging.info("#")
        logging.info("#" + title.upper())
        logging.info("#")
        logging.info(75 * "#")
        logging.info("#")
        
        # capture some information about the progam version, and when and what ran
        logging.info("# Date/Time : " + str(datetime.datetime.now()))
        logging.info("# Python    : " + python_version)
        logging.info("#")
        logging.info("#" + 74 * "-")
        
        # capture the VCF files used for analysis (this should be swapped to the PED file, if used)
        logging.info("# Files")
        logging.info("#" + 74 * "-")
        paths = [self.pedTrio.child.get_path()]
        if self.pedTrio.mother is not None:
            paths.append(self.pedTrio.mother.get_path())
        if self.pedTrio.father is not None:
            paths.append(self.pedTrio.father.get_path())
        for path in paths:
            logging.info("#" + path)
            
        # capture some information about the filters used for screening
        logging.info("#" + 74 * "-")
        logging.info("# Filters")
        logging.info("#" + 74 * "-")
        self.printFileContent(self.filters_path)
        logging.info(75 * "#")
    
    def save_results(self):
        """exports candidate variants and their details (currently matching exome-reporting.pl)
        """
        
        # only add in the header on the first run through
        if self.first_run:
            self.output.write("\t".join(["proband", "alternate_ID", "sex", "chrom", "position", \
                            "gene", "mutation_ID", "transcript", "consequence", "ref/alt_alleles", "MAX_MAF", \
                            "inheritance", "trio_genotype", "mom_aff", "dad_aff", "result"]) + "\n")
        
        if self.pedTrio.father is not None:
            dad_aff = self.pedTrio.father.get_affected_status()
        else:
            dad_aff = "NA"
        if self.pedTrio.mother is not None:
            mom_aff = self.pedTrio.mother.get_affected_status()
        else:
            mom_aff = "NA"
        
        # include an alternate ID for the affected child, if it exists
        if self.ID_mapper is not None:
            alternate_ID = self.ID_mapper[self.pedTrio.child.get_ID()]
        else:
            alternate_ID = 'no_alternate_ID'
        
        reported_some_variants = False
        for candidate in sorted(self.found_variants):
            var = candidate[0]
            result = candidate[1]
            inheritance_type = candidate[2]
            
            # make sure we report the PolyPhen and SIFT scores, if available.
            consequence = var.child.info["CQ"]
            if "PolyPhen" in var.child.info:
                consequence += ",PolyPhen=" + str(var.child.info["PolyPhen"])
            if "SIFT" in var.child.info:
                consequence += ",SIFT=" + str(var.child.info["SIFT"])
            
            transcript = var.child.info["ENST"]
            alleles = var.child.ref_allele + "/" + var.child.alt_allele
            trio_genotype = "%s/%s/%s" % var.get_trio_genotype()
            max_maf = var.child.find_max_allele_frequency(self.tags_dict["MAX_MAF"])
            
            output_line = [self.pedTrio.child.get_ID(), alternate_ID, \
                           self.pedTrio.child.get_gender(), var.get_chrom(), var.get_position(), \
                           var.get_gene(), var.child.get_mutation_id(), transcript, consequence, alleles, \
                           max_maf, inheritance_type, trio_genotype, mom_aff, \
                           dad_aff, result]
            output_line = "\t".join(output_line) + "\n"
            self.output.write(output_line)
            reported_some_variants = True 
        
        # leave a gap between individuals, as per previous reporting system
        if reported_some_variants: 
            self.output.write("\n")
    
    def save_vcf(self):
        """ exports a VCF file for the childs candidate variants.
        """
        
        # figure out what to do with the header
        child_lines = self.vcf_loader.header_lines
        
        child_lines.insert(-1, '##INFO=<ID=ClinicalFilterType,Number=.,Type=String,Description="The type of clinical filter that passed this variant.">\n')
        child_lines.insert(-1, '##INFO=<ID=ClinicalFilterRunDate,Number=.,Type=String,Description="The date on which the clinical filter was run.">\n')
        child_lines.insert(-1, '##INFO=<ID=ClinicalFilterVersion,Number=.,Type=String,Description="The git tag of the clinical filter code.">\n')
        
        child_lines.insert(-1, '##FORMAT=<ID=INHERITANCE_GENOTYPE,Number=.,Type=String,Description="The 012 coded genotypes for a trio (child, mother, father).">\n')
        child_lines.insert(-1, '##FORMAT=<ID=INHERITANCE,Number=.,Type=String,Description="The inheritance of the variant in the trio (biparental, paternal, maternal, deNovo).">\n')
        
        ClinicalFilterRunDate = ",ClinicalFilterRunDate=" + str(datetime.date.today())
        ClinicalFilterVersion = ",ClinicalFilterVersion={0}".format(self.clinicalFilterVersion())
        
        for candidate in sorted(self.found_variants):
            var = candidate[0]
            
            chrom = var.get_chrom()
            position = var.get_position
            
            snp_key = (chrom, position)
            vcf_line = var.child.get_vcf_line()
            
            ClinicalFilterType = ",ClinicalFilterType=" + "XXX"
            vcf_line[7] += ClinicalFilterType + ClinicalFilterRunDate + ClinicalFilterVersion
            
            mother_genotype = var.mother.get_genotype()
            father_genotype = var.father.get_genotype()
            
            if mother_genotype == 0 and father_genotype == 0:
                parental_inheritance = "deNovo"
            elif mother_genotype == 0 and father_genotype != 0:
                parental_inheritance == "paternal"
            elif mother_genotype != 0 and father_genotype == 0:
                parental_inheritance = "maternal"
            else:
                parental_inheritance = "biparental"
            
            vcf_line[8] = ":".join(["INHERITANCE", "INHERITANCE_GENOTYPE"])
            trio_genotype = list(map(str, var.get_trio_genotype()))
            vcf_line[9] = ":".join([parental_inheritance, ",".join(trio_genotype)])
            
            child_lines.append("\t".join(vcf_line) + "\n")
        
        
        # join the list of lines for the VCF file into a single string
        child_lines = "".join(child_lines)
        if platform.python_version_tuple()[0] == "2":
            with gzip.open(self.pedTrio.child.get_ID() + ".vcf.gz", 'wb') as f:
                f.write(child_lines)
        else:
            with gzip.open(self.pedTrio.child.get_ID() + ".vcf.gz", 'wt') as f:
                f.write(child_lines)
        
