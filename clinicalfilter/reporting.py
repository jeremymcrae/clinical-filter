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

import logging
import datetime
import gzip
import importlib
import sys
import os

import clinicalfilter

class Report(object):
    ''' A class to report candidate variants.
    '''
    
    def __init__(self, output_path=None, export_vcf=None, known_genes_date=None):
        ''' initialise the class
        
        Args:
            output_path: path string to list filtered variants in, or None
            export vcf: path string to export VCF files(s), or None
            known_genes_date: date the known gene list was generated, or None
        '''
        
        self.output_path = output_path
        self.export_vcf = export_vcf
        self.known_genes_date = known_genes_date
        
        # clear the tabular output file if it exists
        if self.output_path is not None:
            with open(self.output_path, 'w') as handle:
                handle.write('\t'.join(['proband', 'sex', 'chrom', 'position',
                    'gene', 'mutation_ID', 'transcript', 'consequence',
                    'ref/alt_alleles', 'MAX_MAF', 'inheritance',
                    'trio_genotype', 'mom_aff', 'dad_aff', 'result', 'pp_dnm',
                    'exac_allele_count', 'GQ', 'has_parents', 'cnv_length']) + '\n')
        
        self._log_run_details()
    
    def _log_run_details(self):
        ''' log the python version and run date
        '''
        
        # capture some information about the python version, and run date
        logging.info('# Date/Time : ' + str(datetime.datetime.now()))
        logging.info('# Python    : ' + sys.version.replace('\n', ''))
        logging.info('#')
    
    def export_data(self, variants, family, vcf_header, vcf_provenance):
        ''' export the variants to files (if we have specified paths)
        
        Args:
            variants: list of (variant, check, inheritance) tuples
            family: Family object
        '''
        
        self.family = family
        
        # export the results in tabular format
        if self.output_path is not None:
            self._save_tabular(variants, family)
        
        # export the results in vcf format
        if self.export_vcf is not None:
            lines = self._get_vcf_lines(variants, vcf_header, vcf_provenance)
            path = self._get_vcf_export_path()
            self._write_vcf(path, lines)
    
    def _get_output_line(self, candidate, family):
        ''' gets a tab-separated string for output
        
        Args:
            candidate: (variant, check, inheritance) tuple
            family: Family object for the trio.
        
        Returns:
            tab-separated line in output format
        '''
        
        # get the affected status of the parents
        dad_aff, mom_aff = 'NA', 'NA'
        if family.has_parents():
            dad_aff = family.father.get_affected_status()
            mom_aff = family.mother.get_affected_status()
        
        var = candidate[0]
        
        # make sure we report the PolyPhen and SIFT scores, if available.
        consequence = var.child.info['CQ']
        if 'PolyPhen' in var.child.info:
            consequence += ',PolyPhen={}'.format(var.child.info['PolyPhen'])
        if 'SIFT' in var.child.info:
            consequence += ',SIFT={}'.format(var.child.info['SIFT'])
        if var.child.mnv_code is not None:
            consequence += ',CANDIDATE_MNV={}'.format(var.child.mnv_code)
        
        transcript = 'NA'
        if 'ENST' in var.child.info:
            transcript = var.child.info['ENST']
        
        alleles = '{}/{}'.format(var.child.ref_allele, ','.join(var.child.alt_alleles))
        trio_genotype = '{0}/{1}/{2}'.format(*var.get_trio_genotype())
        trio_genotype = trio_genotype.replace('None', 'NA')
        
        max_af = var.child.find_max_allele_frequency()
        if max_af is None:
            max_af = 'NA'
        max_af = str(max_af)
        
        pp_dnm = 'NA'
        if 'PP_DNM' in var.child.format:
            pp_dnm = var.child.format['PP_DNM']
        
        exac_ac = 'NA'
        if 'AC_Adj' in var.child.info:
            exac_ac = var.child.info['AC_Adj']
        
        cnv_length = 'NA'
        if var.is_cnv():
            start, end = var.get_range()
            cnv_length = str(end - start)
        
        gq = 'NA'
        if 'GQ' in var.child.format:
            gq = str(var.child.format['GQ'])
        
        genes = ','.join(list(set(candidate[3])))
        result = ','.join(candidate[1])
        inh = ','.join(candidate[2])
        
        output_line = [family.child.get_id(), family.child.get_gender(),
            var.get_chrom(), str(var.get_position()), genes,
            var.child.get_mutation_id(), transcript, consequence, alleles,
            max_af, inh, trio_genotype, mom_aff, dad_aff, result, pp_dnm,
            exac_ac, gq, str(family.has_parents()), cnv_length]
        
        return '\t'.join(output_line) + '\n'
    
    def _save_tabular(self, variants, family):
        ''' exports candidate variants and their details
        
        Args:
            variants: list of (variant, check, inheritance) tuples
        '''
        
        # get the affected status of the parents
        dad_aff, mom_aff = 'NA', 'NA'
        if family.has_parents():
            dad_aff = family.father.get_affected_status()
            mom_aff = family.mother.get_affected_status()
        
        with open(self.output_path, 'a') as handle:
            for var in sorted(variants):
                line = self._get_output_line(var, family)
                handle.write(line)
    
    def _get_provenance(self, provenance, member):
        ''' gets the VCF filename, checksum and VCF date for family members
        
        Args:
            provenance: (checksum, path, date) tuple for VCF file
            member: code for member (eg 'proband', 'maternal', 'paternal')
        
        Returns:
            list of lines to add to VCF file
        '''
        
        checksum, path, date = provenance
        sample_id = path.split('.')[0]
        
        sample_id = '##UberVCF_{}_Id={}\n'.format(member, sample_id)
        checksum = '##UberVCF_{}_Checksum={}\n'.format(member, checksum)
        basename = '##UberVCF_{}_Basename={}\n'.format(member, path)
        date = '##UberVCF_{}_Date={}\n'.format(member, date)
        
        return [sample_id, checksum, basename, date]
    
    def _get_vcf_export_path(self):
        ''' get the path for writing a VCF file
        
        Since we optionally define a folder, or path for exporting, we need to
        figure out the path to export a VCF file to.
        
        Returns:
            path to write a vcf file to
        '''
        
        vcf_path = self.export_vcf
        proband_filename = self.family.child.get_id() + ".vcf.gz"
        # check if we have named what looks like a VCF file
        if 'vcf' in vcf_path[-7:] or vcf_path.endswith('gz'):
            # make sure we haven't named a nonexistent folder
            if not os.path.lexists(os.path.dirname(vcf_path)):
                raise ValueError('Cannot find the folder to place the VCF in')
        # if we have named a folder path, add the proband ID for the filename
        elif os.path.isdir(vcf_path):
            vcf_path = os.path.join(vcf_path, proband_filename)
        else:
            raise ValueError('Cannot find the path to export the VCF file')
        
        return vcf_path
    
    def _make_vcf_header(self, header, vcf_provenance):
        ''' start a vcf header using the proband's header, and add extra lines
        
        Args:
            header: list of header lines from the proband's VCF file
            vcf_provenance: list of (checksum, path, date) tuples for family
        
        Returns:
            list of vcf header lines
        '''
        
        # get the final line, then drop it out, so we can insert other lines
        final_header_line = header[-1]
        header = header[:-1]
        
        # define the flags that we add to the info and format fields
        header.append('##INFO=<ID=ClinicalFilterType,Number=.,Type=String,'
            'Description="The type of clinical filter that passed this '
            'variant.">\n')
        header.append('##INFO=<ID=ClinicalFilterGeneInheritance,Number=.,'
            'Type=String,Description="The inheritance mode (Monoallelic, '
            'Biallelic etc) under which the variant was found.">\n')
        header.append('##INFO=<ID=ClinicalFilterReportableHGNC,Number=.,'
            'Type=String,Description="The HGNC symbol which the variant was '
            'identified as being reportable for.">\n')
        header.append('##FORMAT=<ID=INHERITANCE_GENOTYPE,Number=.,Type=String,'
            'Description="The 012 coded genotypes for a trio (child, mother, '
            'father).">\n')
        header.append('##FORMAT=<ID=INHERITANCE,Number=.,Type=String,'
            'Description="The inheritance of the variant in the trio '
            '(biparental, paternal, maternal, deNovo).">\n')
        
        header.append('##ClinicalFilterRunDate={0}\n'.format(datetime.date.today()))
        header.append('##ClinicalFilterVersion={0}\n'.format(clinicalfilter.__version__))
        
        filter_list = ['single_variant', 'compound_het']
        header.append('##ClinicalFilterHistory={0}\n'.format(','.join(filter_list)))
        
        if self.known_genes_date is not None:
            header.append('##ClinicalFilterKnownGenesDate={0}\n'.format(self.known_genes_date))
        
        # add details of the input VCF files used for filtering
        header += self._get_provenance(vcf_provenance[0], 'proband')
        header += self._get_provenance(vcf_provenance[1], 'maternal')
        header += self._get_provenance(vcf_provenance[2], 'paternal')
        
        # add the final header line back in
        header.append(final_header_line)
        
        return header
    
    def _get_parental_inheritance(self, var):
        ''' figures out the parental inheritance for SNVs
        
        Args:
            var: TrioGenotypes object
        
        Returns:
            string for how the variant is inherited eg biparental, deNovo,
            paternal or maternal
        '''
        
        if self.family.has_parents():
            mother_genotype = var.mother.get_genotype()
            father_genotype = var.father.get_genotype()
            
            parental_inheritance = 'biparental'
            if mother_genotype == 0 and father_genotype == 0:
                parental_inheritance = 'deNovo'
            elif mother_genotype == 0 and father_genotype != 0:
                parental_inheritance = 'paternal'
            elif mother_genotype != 0 and father_genotype == 0:
                parental_inheritance = 'maternal'
        else:
            parental_inheritance = 'unknown'
        
        return parental_inheritance
    
    def _get_vcf_lines(self, variants, header, vcf_provenance):
        ''' gets the VCF lines for the proband, including candidate variants.
        
        Args:
            variants: list of (variant, check, inheritance) tuples
            header: list of header lines from the proband's VCF file
            vcf_provenance: list of (checksum, path, date) tuples for family
        
        Returns:
            full list of lines for a VCF file
        '''
        
        vcf_lines = self._make_vcf_header(header, vcf_provenance)
        
        for candidate in sorted(variants):
            var = candidate[0]
            
            vcf_line = var.child.get_vcf_line()
            
            if var.child.mnv_code is not None:
                var.child.add_info_field('CANDIDATE_MNV', var.child.mnv_code)
            
            var.child.add_info_field('ClinicalFilterType', ','.join(candidate[1]))
            var.child.add_info_field('ClinicalFilterGeneInheritance',
                ','.join(candidate[2]))
            var.child.add_info_field('ClinicalFilterReportableHGNC',
                ','.join(list(set(candidate[3]))))
            vcf_line[7] = var.child.get_info_as_string()
            
            parental_inheritance = self._get_parental_inheritance(var)
            
            if 'INHERITANCE' not in vcf_line[8]:
                vcf_line[8] += ':INHERITANCE'
                vcf_line[9] += ':' + parental_inheritance
            
            if not var.is_cnv():
                trio_genotype = '{0},{1},{2}'.format(*var.get_trio_genotype())
                trio_genotype = trio_genotype.replace('None', 'NA')
                vcf_line[8] += ':INHERITANCE_GENOTYPE'
                vcf_line[9] += ':' + trio_genotype
            
            vcf_lines.append('\t'.join(vcf_line) + '\n')
        
        return vcf_lines
    
    def _write_vcf(self, path, lines):
        ''' writes a set of lines to a gzip file
        
        Args:
            path: path to write a file to
            lines: list of lines for a VCF file
        '''
        
        with gzip.open(path, 'w') as handle:
            for x in lines:
                handle.write(x.encode('utf8'))
        
