""" class for reporting found variants via text or VCF files
"""

import logging
import datetime
import gzip
import importlib
import sys
import os

class Report(object):
    """ A class to report candidate variants.
    """
    
    def __init__(self, output_path, export_vcf, ID_mapper, known_genes_date=None):
        """ initialise the class
        
        Args:
            output_path: path string to list filtered variants in, or None
            export vcf: path string to export VCF files(s), or None
            ID_mapper: original_ID - alternate ID dictionary for study probands
            known_genes_date: date the known gene list was generated, or None
        """
        
        self.output_path = output_path
        self.export_vcf = export_vcf
        self.ID_mapper = ID_mapper
        self.known_genes_date = known_genes_date
        
        # clear the tabular output file if it exists
        if self.output_path is not None:
            output = open(self.output_path, "w")
            output.write("\t".join(["proband", "alternate_ID", "sex", \
                "chrom", "position", "gene", "mutation_ID", "transcript", \
                "consequence", "ref/alt_alleles", "MAX_MAF", "inheritance", \
                "trio_genotype", "mom_aff", "dad_aff", "result", "pp_dnm"]) + "\n")
            output.close()
        
        self._log_run_details()
    
    def _clinicalFilterVersion(self):
        """ get the version (git tag) from clinicalfilter.version.version()
        
        Returns:
            string for the git tag, or "XXX"
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
    
    def _log_run_details(self):
        """ log the python version and run date
        """
        
        # capture the program title
        logging.info("# Clinical filtering analysis")
        
        # capture some information about the python version, and run date
        logging.info("# Date/Time : " + str(datetime.datetime.now()))
        logging.info("# Python    : " + sys.version.replace("\n", ""))
        logging.info("#")
    
    def export_data(self, variants, family, vcf_header, vcf_provenance):
        """ export the variants to files (if we have specified paths)
        
        Args:
            variants: list of (variant, check, inheritance) tuples
            family: Family object
        """
        
        self.family = family
        
        # export the results in tabular format
        if self.output_path is not None:
            self.output = open(self.output_path, "a")
            self._save_tabular(variants)
            self.output.close()
        
        # export the results in vcf format
        if self.export_vcf is not None:
            vcf_lines = self._get_vcf_lines(variants, vcf_header, vcf_provenance)
            vcf_path = self._get_vcf_export_path()
            self._write_vcf(vcf_path, vcf_lines)
    
    def _get_output_line(self, candidate, dad_aff, mom_aff, alternate_ID):
        """ gets a tab-separated string for output
        
        Args:
            candidate: (variant, check, inheritance) tuple
            dad_aff: affected status of father (or "NA")
            mom_aff: affected status of mother (or "NA")
            alternate_ID: alternate ID for proband (or "NA")
        
        Returns:
            tab-separated line in output format
        """
        
        var = candidate[0]
        
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
        trio_genotype = "{0}/{1}/{2}".format(*var.get_trio_genotype())
        
        max_maf = var.child.find_max_allele_frequency()
        if max_maf is None:
            max_maf = "NA"
        max_maf = str(max_maf)
        
        pp_dnm = "NA"
        if "PP_DNM" in var.child.format:
            pp_dnm = var.child.format["PP_DNM"]
        
        genes = ",".join(list(set(candidate[3])))
        result = ",".join(candidate[1])
        inh = ",".join(candidate[2])
        
        output_line = [self.family.child.get_id(), alternate_ID, \
            self.family.child.get_gender(), var.get_chrom(), \
            str(var.get_position()), genes, var.child.get_mutation_id(), \
            transcript, consequence, alleles, max_maf, inh, \
            trio_genotype, mom_aff, dad_aff, result, pp_dnm]
        
        output_line = "\t".join(output_line) + "\n"
        
        return output_line
    
    def _save_tabular(self, variants):
        """ exports candidate variants and their details
        
        Args:
            variants: list of (variant, check, inheritance) tuples
        """
        
        # get the affected status of the parents
        dad_aff = "NA"
        mom_aff = "NA"
        if self.family.has_parents():
            dad_aff = self.family.father.get_affected_status()
            mom_aff = self.family.mother.get_affected_status()
        
        # include an alternate ID for the affected child, if it exists
        alt_id = 'no_alternate_ID'
        if self.ID_mapper is not None:
            alt_id = self.ID_mapper[self.family.child.get_id()]
        
        for var in sorted(variants):
            output_line = self._get_output_line(var, dad_aff, mom_aff, alt_id)
            self.output.write(output_line)
        
        # leave a gap between individuals, as per previous reporting system
        if len(variants) > 0:
            self.output.write("\n")
    
    def _get_provenance(self, provenance, member):
        """ gets the VCF filename, checksum and VCF date for family members
        
        Args:
            provenance: (checksum, path, date) tuple for VCF file
            member: code for member (eg "proband", "maternal", "paternal")
        
        Returns:
            list of lines to add to VCF file
        """
        
        ID = "##UberVCF_" + member + "_Id=" + provenance[1].split(".")[0] + "\n"
        checksum = "##UberVCF_" + member + "_Checksum=" + provenance[0] + "\n"
        basename = "##UberVCF_" + member + "_Basename=" + provenance[1] + "\n"
        date = "##UberVCF_" + member + "_Date=" + provenance[2] + "\n"
        
        return [ID, checksum, basename, date]
    
    def _get_vcf_export_path(self):
        """ get the path for writing a VCF file
        
        Since we optionally define a folder, or path for exporting, we need to
        figure out the path to export a VCF file to.
        
        Returns:
            path to write a vcf file to
        """
        
        vcf_path = self.export_vcf
        proband_filename = self.family.child.get_id() + ".vcf.gz"
        # check if we have named what looks like a VCF file
        if "vcf" in vcf_path[-7:] or vcf_path.endswith("gz"):
            # make sure we haven't named a nonexistent folder
            if not os.path.lexists(os.path.dirname(vcf_path)):
                raise ValueError("Cannot find the folder to place the VCF in")
        # if we have named a folder path, add the proband ID for the filename
        elif os.path.isdir(vcf_path):
            vcf_path = os.path.join(vcf_path, proband_filename)
        else:
            raise ValueError("Cannot find the path to export the VCF file")
        
        return vcf_path
    
    def _make_vcf_header(self, header, vcf_provenance):
        """ start a vcf header using the proband's header, and add extra lines
        
        Args:
            header: list of header lines from the proband's VCF file
            vcf_provenance: list of (checksum, path, date) tuples for family
        
        Returns:
            list of vcf header lines
        """
        
        # get the final line, then drop it out, so we can insert other lines
        final_header_line = header[-1]
        header = header[:-1]
        
        # define the flags that we add to the info and format fields
        header.append('##INFO=<ID=ClinicalFilterType,Number=.,Type=String,Description="The type of clinical filter that passed this variant.">\n')
        header.append('##INFO=<ID=ClinicalFilterGeneInheritance,Number=.,Type=String,Description="The inheritance mode (Monoallelic, Biallelic etc) under which the variant was found.">\n')
        header.append('##INFO=<ID=ClinicalFilterReportableHGNC,Number=.,Type=String,Description="The HGNC symbol which the variant was identified as being reportable for.">\n')
        header.append('##FORMAT=<ID=INHERITANCE_GENOTYPE,Number=.,Type=String,Description="The 012 coded genotypes for a trio (child, mother, father).">\n')
        header.append('##FORMAT=<ID=INHERITANCE,Number=.,Type=String,Description="The inheritance of the variant in the trio (biparental, paternal, maternal, deNovo).">\n')
        
        header.append("##ClinicalFilterRunDate={0}\n".format(datetime.date.today()))
        header.append("##ClinicalFilterVersion={0}\n".format(self._clinicalFilterVersion()))
        
        filter_list = ["single_variant", "compound_het"]
        header.append("##ClinicalFilterHistory={0}\n".format(",".join(filter_list)))
        
        if self.known_genes_date is not None:
            header.append("##ClinicalFilterKnownGenesDate={0}\n".format(self.known_genes_date))
        
        # add details of the input VCF files used for filtering
        header += self._get_provenance(vcf_provenance[0], "proband")
        header += self._get_provenance(vcf_provenance[1], "maternal")
        header += self._get_provenance(vcf_provenance[2], "paternal")
        
        # add the final header line back in
        header.append(final_header_line)
        
        return header
    
    def _get_parental_inheritance(self, var):
        """ figures out the parental inheritance for SNVs
        
        Args:
            var: TrioGenotypes object
        
        Returns:
            string for how the variant is inherited eg biparental, deNovo,
            paternal or maternal
        """
        
        if self.family.has_parents():
            mother_genotype = var.mother.get_genotype()
            father_genotype = var.father.get_genotype()
            
            parental_inheritance = "biparental"
            if mother_genotype == 0 and father_genotype == 0:
                parental_inheritance = "deNovo"
            elif mother_genotype == 0 and father_genotype != 0:
                parental_inheritance = "paternal"
            elif mother_genotype != 0 and father_genotype == 0:
                parental_inheritance = "maternal"
        else:
            parental_inheritance = "unknown"
        
        return parental_inheritance
    
    def _get_vcf_lines(self, variants, header, vcf_provenance):
        """ gets the VCF lines for the proband, including candidate variants.
        
        Args:
            variants: list of (variant, check, inheritance) tuples
            header: list of header lines from the proband's VCF file
            vcf_provenance: list of (checksum, path, date) tuples for family
        
        Returns:
            full list of lines for a VCF file
        """
        
        vcf_lines = self._make_vcf_header(header, vcf_provenance)
        
        for candidate in sorted(variants):
            var = candidate[0]
            
            vcf_line = var.child.get_vcf_line()
            
            filter_type = ";ClinicalFilterType=" + ",".join(candidate[1])
            gene_inheritance = ";ClinicalFilterGeneInheritance=" + ",".join(candidate[2])
            reportable_gene = ";ClinicalFilterReportableHGNC={0}".format(",".join(list(set(candidate[3]))))
            vcf_line[7] += gene_inheritance + filter_type + reportable_gene
            
            parental_inheritance = self._get_parental_inheritance(var)
            
            if "INHERITANCE" not in vcf_line[8]:
                vcf_line[8] += ":INHERITANCE"
                vcf_line[9] += ":" + parental_inheritance
            
            if not var.is_cnv():
                trio_genotype = "{0},{1},{2}".format(*var.get_trio_genotype())
                vcf_line[8] += ":INHERITANCE_GENOTYPE"
                vcf_line[9] += ":" + trio_genotype
            
            vcf_lines.append("\t".join(vcf_line) + "\n")
        
        return vcf_lines
    
    def _write_vcf(self, path, vcf_lines):
        """ writes a set of lines to a gzip file
        
        Args:
            path: path to write a file to
            vcf_lines: list of lines for a VCF file
        """
        
        # join the list of lines for the VCF file into a single string
        vcf_lines = "".join(vcf_lines)
        
        # get a gzip file handle (needs to be python version specific)
        if sys.version_info[0] == 2:
            f = gzip.open(path, 'wb')
        else:
            f = gzip.open(path, 'wt')
        
        f.write(vcf_lines)
        f.close()
        
