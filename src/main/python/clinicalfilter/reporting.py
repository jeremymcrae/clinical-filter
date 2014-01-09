""" class for reporting found variants via text or VCF files
"""


import os
import sys
import logging
import platform
import datetime
import gzip

import user


class report(object):
    """A class to report candidate variants.
    """
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
        user.printFileContent(self.filters_path)
        logging.info(75 * "#")
    
    def save_results(self):
        """exports candidate variants and their details (currently matching exome-reporting.pl)
        """
        results_dict = self.found_variants
        
        # only add in the header on the first run through
        if self.first_run:
            self.output.write("\t".join(["proband", "alternate_ID", "sex", "chrom", "position", \
                            "gene", "transcript", "consequence", "ref/alt_alleles", "MAX_MAF", \
                            "inheritance", "trio_genotype", "mom_aff", "dad_aff", "result"]) + "\n")
        
        reported_some_variants = False
        for gene in sorted(results_dict):
            for position in sorted(results_dict[gene]):
                snp = results_dict[gene][position]
                
                # make sure we report the PolyPhen and SIFT scores, if available.
                consequence = snp["consequence"]
                if snp["PolyPhen"] is not None:
                    consequence += ",PolyPhen=" + str(snp["PolyPhen"])
                if snp["SIFT"] is not None:
                    consequence += ",SIFT=" + str(snp["SIFT"])
                snp["consequence"] = consequence
                
                alleles = snp["REF"] + "/" + snp["ALT"]
                
                output_line = [snp["person_ID"], snp["alternate_ID"], \
                               self.pedTrio.child.get_gender(), snp["CHROM"], snp["POS"], \
                               snp["gene"], snp['transcript'], snp["consequence"], alleles, \
                               snp["MAX_MAF"], snp["inh"], snp["trio_genotype"], snp["mom_aff"], \
                               snp["dad_aff"], snp['result']]
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
        child_lines = self.child_vcf["header_lines"]
        
        child_lines.insert(-1, '##INFO=<ID=ClinicalFilterType,Number=.,Type=String,Description="The type of clinical filter that passed this variant.">\n')
        child_lines.insert(-1, '##INFO=<ID=ClinicalFilterRunDate,Number=.,Type=String,Description="The date on which the clinical filter was run.">\n')
        child_lines.insert(-1, '##INFO=<ID=ClinicalFilterVersion,Number=.,Type=String,Description="The git tag of the clinical filter code.">\n')
        
        child_lines.insert(-1, '##FORMAT=<ID=INHERITANCE_GENOTYPE,Number=.,Type=String,Description="The 012 coded genotypes for a trio (child, mother, father).">\n')
        child_lines.insert(-1, '##FORMAT=<ID=INHERITANCE,Number=.,Type=String,Description="The inheritance of the variant in the trio (biparental, paternal, maternal, deNovo).">\n')
        
        ClinicalFilterRunDate = ",ClinicalFilterRunDate=" + str(datetime.date.today())
        ClinicalFilterVersion = ",ClinicalFilterVersion=XXX"
        
        results_dict = self.found_variants
        for gene in sorted(results_dict):
            for position in sorted(results_dict[gene]):
                snp = results_dict[gene][position]
                chrom = snp["CHROM"]
                position = snp["POS"]
                
                snp_key = (chrom, position)
                vcf_line = self.child_vcf["data"][snp_key]["vcf_line"]
                
                ClinicalFilterType = ",ClinicalFilterType=" + "XXX"
                vcf_line[7] += ClinicalFilterType + ClinicalFilterRunDate + ClinicalFilterVersion
                
                mother_genotype = snp["trio_genotype"].split("/")[1]
                father_genotype = snp["trio_genotype"].split("/")[2]
                
                if mother_genotype == 0 and father_genotype == 0:
                    parental_inheritance = "deNovo"
                elif mother_genotype == 0 and father_genotype != 0:
                    parental_inheritance == "paternal"
                elif mother_genotype != 0 and father_genotype == 0:
                    parental_inheritance = "maternal"
                else:
                    parental_inheritance = "biparental"
                
                vcf_line[8] = ":".join(["INHERITANCE", "INHERITANCE_GENOTYPE"])
                vcf_line[9] = ":".join([parental_inheritance, snp["trio_genotype"].replace("/", ",")])
                
                child_lines.append("\t".join(vcf_line) + "\n")
        
        # join the list of lines for the VCF file into a single string
        child_lines = "".join(child_lines)
        with gzip.open(self.pedTrio.child.get_ID() + ".vcf.gz", 'wb') as f:
            f.write(child_lines)
