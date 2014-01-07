'''
Methods to manipulate vcf file
3 Dec 2011 version 0.1
sa9@sanger.ac.uk
'''

import os
import io
import platform
import sys
import gzip
import logging

import variant


def getVcfFile(path):
    """Gets a VCF file handle while allowing for gzipped and text VCF formats.

    Args:
        path: path to VCF file.

    Returns:
        A file handle for VCF file.

    Raises:
        ValueError: An error when the VCF file path is not specified correctly.
    """

    if not os.path.exists(path):
        raise IOError("VCF file not found at: " + path)

    if path.endswith(".gz"):
        # python2 gzip opens in text, but same mode in python3 opens as bytes,
        # avoid with platform specific code
        if platform.python_version_tuple()[0] == "2":
            f = gzip.open(path, 'r')
        else:
            f = gzip.open(path, "rt")
    elif path.endswith(".vcf") or path.endswith(".txt"):
        f = io.open(path,'r', encoding="latin_1")
    else:
        print('Unable to open files with extension "%s". Accepted extensions are ".vcf" or ".txt"' % path.split(".")[-1])
        print('The ".txt" is tab-separated file and it can be compressed in gz format.')
        sys.exit(0)
    return f

def get_vcf_header(f):
    '''Get the header lines from a VCF file.

    Args:
        f: file object to read from
    
    Returns:
        a list of lines that start with "#", which are the header lines.
    '''
    header = []
    
    for line in f:
        if not line.startswith("#"):
            break
        header.append(line)
    
    return header

def getHeaderLabels(header):
    '''Get the header labels from a VCF file.

    Read VCF header and get all fields arrangd into three groups
     - The main columns are from chr to filter (columns 0-6)
     - The info tag fields (column 7)
     - The format tag fields (column 8 and their values are in 9)
    These labels are used to create empty data dictionary filled with "NA" to be exported as TSV 
    object (or writen to a file)
    '''
    labels_dict = {
        "m":[], # main
        "i":[], # info
        "f":[]  # format
    }
    for line in header:
        # process the FORMAT and INFO ID labels in the VCF header lines
        if line.startswith("##INFO="):
            ID = line.split("ID=")[1].split(",")[0]
            labels_dict["i"].append(ID)
        elif line.startswith("##FORMAT="):
            ID = line.split("ID=")[1].split(",")[0]
            labels_dict["f"].append(ID)
        # process the VCF column headers
        elif line.startswith("#CHROM"):
                line = line.strip().lstrip("#").split("\t")
                for i, item in enumerate(line):
                    if i <= 6:
                        labels_dict["m"].append(item)

    # count the occurences of each unique label in the dictionary
    counts = {}
    for group in ["m", "i", "f"]:
        for label in labels_dict[group]:
            if label not in counts:
                counts[label] = 0
            counts[label] += 1
    
    # create one list and re-name any duplicate labels with its parent letter (e.g. DP_f and DP_i )
    labels = [] # final
    for group in ["m", "i", "f"]:
        for label in labels_dict[group]:
            if counts[label] > 1:
                label = label + "_" + group
            labels.append(label)
    
    return labels

class MatchCNVs(object):
    """ class to find if a parents CNV matches any of a childs CNVs
    """
    
    def __init__(self, child_variants):
        """ initiate the class with the childs variants
        """
        
        self.child_variants = child_variants
        
        # figure out which of the childs variants are CNVs
        self.cnvs = []
        for key in self.child_variants:
            if len(key) == 3: # ignore SNVs, which are only (chrom, position)
                self.cnvs.append(key)
    
    def has_match(self, var):
        """ checks if a CNV has any matches amongst the childs CNVs
        
        Args:
            var: CNV object
        
        Returns:
            returns true if any of the childs CNVs are good matches
        """
        
        # TODO: swap to code that compares CNV sizes in line below
        # if self.any_overlap(var) and self.similar_size(var):
        if self.any_overlap(var):
            print(var.get_key(), self.overlap)
            return True
        else:
            return False
        
    def any_overlap(self, var):
        """ checks if any of the childs CNVs overlap the current CNV
        
        Args:
            var: CNV object
        
        Returns:
            returns true if any of the childs CNVs overlap
        """
        
        var_key = var.get_key() # CNVs are (chrom, start, end)
        var_chrom = var_key[0]
        var_start = int(var_key[1])
        var_end = int(var_key[2])
        
        self.overlap = []
        for key in self.cnvs:
            
            child_chrom = key[0]
            child_start = int(key[1])
            child_end = int(key[2])
            
            if var_chrom != child_chrom:
                continue
            
            # check if the childs CNVs end points lie within the CNV that we
            # are examining, or if the childs CNV surrounds the examined CNV
            if var_end >= child_start >= var_start or \
               var_end >= child_end >= var_start or \
               child_end >= var_start >= child_start:
                self.overlap.append(key)
        
        if len(self.overlap) > 0:
            return True
        else:
            return False
        
    def similar_size(self, var):
        """ checks if the current CNV matches the overlapping child CNVs size
        
        Args:
            var: CNV object
        
        Returns:
            returns true if any of the overlapping CNVs have similar sizes
        """
        (min_size, max_size) = var.calculate_cnv_size_tolerance()
        
        for overlap in self.overlap:
            var_chrom = overlap[0]
            var_start = int(overlap[1])
            var_end = int(overlap[2])
            var_size = var_end - var_start
            
            if max_size > var_size > min_size:
                self.overlap = overlap
                return True
        
        return False
    
    def get_overlap_key(self):
        """ returns the tuple for an overlapping CNV
        """
        
        return self.overlap


def vcf2tsv(path, gender, filters=None, child_variants=None):
    ''' Convert VCF to TSV format. Use for single sample VCF file.

    Obtains the VCF data for a single sample. This function optionally filters the lines of the 
    VCF file that pass defined criteria, in order to reduce memory usage.

    Args:
        path: path to VCF file for the sample.
        gender: gender of the individual
        filters: optional filters, as loaded by parseFilters in user.py
        child_variants: list of tuple keys for variants identified in affected child, in order to
                        get the same variants from the parents.

    Returns:
        A dictionary containing header data, and variant data for each variant.
    '''
    
    if child_variants is not None:
        cnv_matcher = MatchCNVs(child_variants)
    
    try:
        f = getVcfFile(path)
    except IOError as e:
        raise e
    
    header = get_vcf_header(f)
    
    labels = getHeaderLabels(header)
    vcf = {"header": [], "data": {}}
    vcf["header"] = labels
    vcf["header_lines"] = header
    
    for line in f:
        if line.startswith("#"):
            continue
        
        line = line.strip().split("\t")
        
        chrom = line[0]
        position = line[1]
        snp_id = line[2]
        ref_allele = line[3]
        alt_allele = line[4]
        quality = line[5]
        filter_value = line[6]
        info_values = line[7]
        format_keys = line[8]
        sample_values = line[9]
        
        if alt_allele == "<DUP>" or alt_allele == "<DEL>":
            var = variant.CNV(chrom, position, snp_id, ref_allele, alt_allele, quality, filter_value)
            var.add_info(info_values)
        else:
            var = variant.SNV(chrom, position, snp_id, ref_allele, alt_allele, quality, filter_value)
        
        include_variant = False
        if child_variants is not None:
            key = var.get_key()
            if key in child_variants:
                include_variant = True
            elif alt_allele == "<DUP>" or alt_allele == "<DEL>":
                if cnv_matcher.has_match(var):
                    include_variant = True
        elif filters is not None:
            var.add_info(info_values)
            if var.passes_filters(filters):
                include_variant = True
        else:
            include_variant = True
        
        if include_variant:
            var.add_info(info_values)
            var.add_format(format_keys, sample_values)
            var.add_vcf_line(line)
            var.set_gender(gender)
            var.set_genotype()
            
            key = var.get_key()
            vcf["data"][key] = var
        
        # data_holder = dict.fromkeys(labels, None)
        # for i in range(0, 7): #i, item in enumerate(line):
        #     data_holder[labels[i]] = line[i]
        # 
        # # figure out whether to include the variant or not
        # include_variant = False
        # if child_variants is not None:
        #     key = (line[0], line[1])
        #     if key in child_variants:
        #         include_variant = True
        # elif isPassUserFilters is not None:
        #     parseINFO(data_holder, line[7])
        #     if isPassUserFilters(data_holder):
        #         include_variant = True
        # else:
        #     include_variant = True
        # 
        # if include_variant:
        #     key = (line[0], line[1])
        #     parseINFO(data_holder, line[7])
        #     parseFORMAT(data_holder, line[8:])
        #     vcf["data"][key] = data_holder
        #     vcf["data"][key]["vcf_line"] = line
    
    logging.debug(path)
    for key in sorted(vcf["data"]):
        logging.debug(vcf["data"][key])
    
    return vcf

