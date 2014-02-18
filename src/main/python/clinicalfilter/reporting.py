""" class for reporting found variants via text or VCF files
"""

import logging
import platform
import datetime
import gzip
import importlib
import sys
import os

class Report(object):
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
        f = open(path, "r")
        for line in f:
            if not line.startswith("#"):
                logging.info("#" + line.strip())
        f.close()
    
    def printReportInputInfo(self, title):
        """report the parameters used to run the script, such as python version, script version, 
        filters and weights, as well as which columns were chosen for output. 
        """
        
        python_version = "%s %s %s" % (platform.python_version(), platform.python_build(), platform.python_compiler())
        
        # capture the program title
        logging.info("#" + title)
        
        # capture some information about the progam version, and when and what ran
        logging.info("# Date/Time : " + str(datetime.datetime.now()))
        logging.info("# Python    : " + python_version)
        logging.info("#")
            
        # capture some information about the filters used for screening
        logging.info("#" + 74 * "-")
        logging.info("# Filters")
        logging.info("#" + 74 * "-")
        self.printFileContent(self.filters_path)
        logging.info(75 * "#")
    
    def save_results(self, found_vars):
        """exports candidate variants and their details
        """
        
        # only add in the header on the first run through
        if self.first_run:
            self.output.write("\t".join(["proband", "alternate_ID", "sex", \
                "chrom", "position", "gene", "mutation_ID", "transcript", \
                "consequence", "ref/alt_alleles", "MAX_MAF", "inheritance", \
                "trio_genotype", "mom_aff", "dad_aff", "result"]) + "\n")
        
        if self.family.has_parents():
            dad_aff = self.family.father.get_affected_status()
            mom_aff = self.family.mother.get_affected_status()
        else:
            dad_aff = "NA"
            mom_aff = "NA"
        
        # include an alternate ID for the affected child, if it exists
        if self.ID_mapper is not None:
            alternate_ID = self.ID_mapper[self.family.child.get_ID()]
        else:
            alternate_ID = 'no_alternate_ID'
        
        reported_some_variants = False
        for candidate in sorted(found_vars):
            var = candidate[0]
            filter_type = candidate[1]
            inheritance_type = candidate[2]
            
            # make sure we report the PolyPhen and SIFT scores, if available.
            consequence = var.child.info["CQ"]
            if "PolyPhen" in var.child.info:
                consequence += ",PolyPhen=" + str(var.child.info["PolyPhen"])
            if "SIFT" in var.child.info:
                consequence += ",SIFT=" + str(var.child.info["SIFT"])
            
            transcript = "NA"
            if "ENST" in var.child.info:
                transcript = var.child.info["ENST"]
            alleles = var.child.ref_allele + "/" + var.child.alt_allele
            trio_genotype = "%s/%s/%s" % var.get_trio_genotype()
            max_maf = var.child.find_max_allele_frequency(self.tags_dict["MAX_MAF"])
            
            output_line = [self.family.child.get_ID(), alternate_ID, \
                self.family.child.get_gender(), var.get_chrom(), \
                var.get_position(), var.get_gene(),var.child.get_mutation_id(),\
                transcript, consequence, alleles, max_maf, inheritance_type, \
                trio_genotype, mom_aff, dad_aff, filter_type]
            output_line = "\t".join(output_line) + "\n"
            self.output.write(output_line)
            reported_some_variants = True 
        
        # leave a gap between individuals, as per previous reporting system
        if reported_some_variants: 
            self.output.write("\n")
    
    def include_vcf_provenance(self, provenance, member):
        """ adds the original VCF filename and checksums
        
        Args:
            provenance: tuple of checksum, filename and date for a member of a trio_genotype
            member: code for member (eg "proband", "maternal", "paternal")
        
        Returns:
            list of lines to add to VCF file
        """
        
        ID = "##UberVCF_" + member + "_Id=" + provenance[1].split(".")[0] + "\n"
        checksum = "##UberVCF_" + member + "_Checksum=" + provenance[0] + "\n"
        basename = "##UberVCF_" + member + "_Basename=" + provenance[1] + "\n"
        date = "##UberVCF_" + member + "_Date=" + provenance[2] + "\n"
        
        return [ID, checksum, basename, date]
    
    def get_vcf_export_path(self):
        """ get the path for writing a VCF file
        """
        
        vcf_path = self.export_vcf
        proband_filename = self.family.child.get_ID() + ".vcf.gz"
        # check if we have named what looks like a VCF file
        if "vcf" in vcf_path[-7:] or vcf_path.endswith("gz"):
            # make sure we haven't named a nonexistent folder
            if not os.path.lexists(os.path.dirname(vcf_path)):
                raise ValueError("Cannot find the folder to export the VCF file")
        # if we have named a folder path, add the proband ID for the filename
        elif os.path.isdir(vcf_path):
            vcf_path = os.path.join(vcf_path, proband_filename)
        else:
            raise ValueError("Cannot find the path to export the VCF file")
        
        return vcf_path
    
    def save_vcf(self, found_vars):
        """ exports a VCF file for the childs candidate variants.
        """
        
        # figure out what to do with the header
        child_lines = self.vcf_loader.child_header
        
        child_lines.insert(-1, '##INFO=<ID=ClinicalFilterType,Number=.,Type=String,Description="The type of clinical filter that passed this variant.">\n')
        child_lines.insert(-1, '##INFO=<ID=ClinicalFilterGeneInheritance,Number=.,Type=String,Description="The inheritance mode (Monoallelic, Biallelic etc) under which the variant was found.">\n')
        # child_lines.insert(-1, '##INFO=<ID=ClinicalFilterRunDate,Number=.,Type=String,Description="The date on which the clinical filter was run.">\n')
        # child_lines.insert(-1, '##INFO=<ID=ClinicalFilterVersion,Number=.,Type=String,Description="The git tag of the clinical filter code.">\n')
        
        child_lines.insert(-1, '##FORMAT=<ID=INHERITANCE_GENOTYPE,Number=.,Type=String,Description="The 012 coded genotypes for a trio (child, mother, father).">\n')
        child_lines.insert(-1, '##FORMAT=<ID=INHERITANCE,Number=.,Type=String,Description="The inheritance of the variant in the trio (biparental, paternal, maternal, deNovo).">\n')
        
        child_lines.insert(-1, "##ClinicalFilterRunDate={0}\n".format(datetime.date.today()))
        child_lines.insert(-1, "##ClinicalFilterVersion={0}\n".format(self.clinicalFilterVersion()))
        
        filter_list = ["single_variant", "compound_het"]
        child_lines.insert(-1, "##ClinicalFilterHistory={0}\n".format(",".join(filter_list)))
        
        if hasattr(self, "known_genes_date"): 
            child_lines.insert(-1, "##ClinicalFilterKnownGenesDate={0}\n".format(self.known_genes_date))
        
        child_lines = child_lines[:-1] + self.include_vcf_provenance(self.vcf_provenance[0], "proband") + child_lines[-1:]
        child_lines = child_lines[:-1] + self.include_vcf_provenance(self.vcf_provenance[1], "maternal") + child_lines[-1:]
        child_lines = child_lines[:-1] + self.include_vcf_provenance(self.vcf_provenance[2], "paternal") + child_lines[-1:]
        
        var_lines = []
        filter_strings = set([])
        for candidate in sorted(found_vars):
            var = candidate[0]
            filter_type = candidate[1]
            gene_inheritance = candidate[2]
            
            chrom = var.get_chrom()
            position = var.get_position
            
            vcf_line = var.child.get_vcf_line()
            
            ClinicalFilterType = ";ClinicalFilterType=" + filter_type
            ClinicalFilterGeneInheritance = ";ClinicalFilterGeneInheritance=" + gene_inheritance
            vcf_line[7] += ClinicalFilterGeneInheritance + ClinicalFilterType
            
            filter_strings.add(filter_type + "," + gene_inheritance)
            
            if self.family.has_parents():
                mother_genotype = var.mother.get_genotype()
                father_genotype = var.father.get_genotype()
                
                parental_inheritance = "biparental"
                if mother_genotype == 0 and father_genotype == 0:
                    parental_inheritance = "deNovo"
                elif mother_genotype == 0 and father_genotype != 0:
                    parental_inheritance == "paternal"
                elif mother_genotype != 0 and father_genotype == 0:
                    parental_inheritance = "maternal"
            else:
                parental_inheritance = "unknown"
            
            if "INHERITANCE" not in vcf_line[8]:
                vcf_line[8].append("INHERITANCE")
                vcf_line[9].append(parental_inheritance)
            
            if not var.child.is_cnv():
                trio_genotype = "%s,%s,%s" % var.get_trio_genotype()
                vcf_line[8].append("INHERITANCE_GENOTYPE")
                vcf_line[9].append(trio_genotype)
            
            vcf_line[8] = ":".join(vcf_line[8])
            vcf_line[9] = ":".join(vcf_line[9])
            
            var_lines.append("\t".join(vcf_line) + "\n")
        
        
        child_lines += var_lines
        
        vcf_path = self.get_vcf_export_path()
        # join the list of lines for the VCF file into a single string
        child_lines = "".join(child_lines)
        if platform.python_version_tuple()[0] == "2":
            with gzip.open(vcf_path, 'wb') as f:
                f.write(child_lines)
        else:
            with gzip.open(vcf_path, 'wt') as f:
                f.write(child_lines)
        
