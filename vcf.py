'''
Methods to manipulate vcf file
3 Dec 2011 version 0.1
sa9@sanger.ac.uk
'''

import os
import sys
import gzip
import logging


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
        f =  gzip.open(path,'r')
    elif path.endswith(".vcf") or path.endswith(".txt"):
        f = open(path,'r')
    else:
        print 'Unable to open files with extension "%s". Accepted extensions are ".vcf" or ".txt"' % path.split(".")[-1]
        print 'The ".txt" is tab-separated file and it can be compressed in gz format.'
        sys.exit(0)
    return f

# def fliePath(f):
#     """ I don't think this function is used anywhere, possibly remove.
#     """
#     return f.name

def getKey(line, type="chr_pos"):
    """Creates a unique key for the variants being investigated.

    Args:
        line: list of values from a split tab-separated VCF line
        type: the kind of unique ID to create from the position and allele codes.

    Returns:
        set pf values defineing the unique identifier for the variant based on the variant position
        and allele codes. For example, this variant:

        chrom: "1"
        pos: "14907"
        reference allele: "A"
        alternate allele: "G"

        returns ("1", "14907", "A", "G")
    """
    if type == "pos":
        key = (line[1])
    elif type == "chr_pos":
        key = (line[0], line[1])
    elif type == "chr_pos_ref_alt":
        key = (line[0], line[1], line[3], line[4]) # chr,pos, ref, alt
    
    #if line[0] == "X":
    #    print key
    return key


def translateGT(GT):
    """Maps genotypes from two character format to single character.

    Args:
        GT: genotype in two character format. eg '0/0'

    Returns:
        Genotype in single character format. eg '0'
    """
    # This function might run quicker as the following. Possibly could speed up further by passing 
    # reference to a dictionary, rather than creating new each time.
    genotype_dict = {'00': 0, '01': 1, '10': 1, '12': 1, '21': 1, '02': 1, '20': 1, '11': 2, '22': 2}
    return genotype_dict[GT[0] + GT[2]]

def parseINFO(info):
    """Parses the INFO column from VCF files.

    Args:
        # data_holder: dictionary entry for a VCF line
        info: INFO text from a line in a VCF file

    Returns:
        A dictionary including INFO items, if the item is a key-value pair, it uses that for the 
        entry, otherwise uses the key and value of '1' for the entry.
    """
    data_holder = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=")
        else:
            k, v = item, "1"
        data_holder[k] = v
    return data_holder

def parseFORMAT(data_holder, formats):
    """Parses the FORMAT column from VCF files.

    Args:
        data_holder: dictionary entry for a VCF line
        formats: FORMAT text from a line in a VCF file

    Returns:
        A dictionary including FORMAT items, which uses each FORMAT label as a key, and the 
         value is the corresponding value in following column. For example:

        {'DP': '10', 'FP': '40', etc}
    """
    tag_labels = formats[0].split(":") # the first item in formats are the tags DP:FP:ETC
    tag_values = formats[1].split(":") # the second item are the corresponding values.
    
    for i, value in enumerate(tag_values):
        data_holder[tag_labels[i]] = value
    return data_holder

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

def vcf2tsv(path, isPassUserFilters=None, child_variants=None):
    ''' Convert VCF to TSV format. Use for single sample VCF file.

    Obtains the VCF data for a single sample. This function optionally filters the lines of the 
    VCF file that pass defined criteria, in order to reduce memory usage.

    Args:
        path: path to VCF file for the sample.
        isPassUserFilters: optional filters, as loaded by parseFilters in user.py
        child_variants: list of tuple keys for variants identified in affected child, in order to
                        get the same variants from the parents.

    Returns:
        A dictionary containing header data, and variant data for each variant.
    '''
    
    try:
        f = getVcfFile(path)
    except IOError as e:
        raise e
    
    header = get_vcf_header(f)
    
    labels = getHeaderLabels(header)
    vcf = {"header": [], "data": {}}
    vcf["header"] = labels
    vcf["header_lines"] = header
    
    # total_variant_count = 0
    # passing_filters = 0
    for line in f:
        if not line.startswith("#"):
        #if line.startswith("X"):  #uncomment this line if you want to test one chromosome
            line = line.strip().split("\t")
            # total_variant_count += 1
            
            # start the
            data_holder = dict.fromkeys(labels, None)
            data_holder.update(parseINFO(line[7]))
            for i in range(0, 7): #i, item in enumerate(line):
                data_holder[labels[i]] = line[i]
            
            # if the isPassUserFilters function is passed into this function, then test the current
            # record for whether it passes the filtering criteria such as low minor allele 
            # frequency etc
            if child_variants is not None:
                key = getKey(line, "chr_pos")
                if key in child_variants:
                    data_holder = parseFORMAT(data_holder, line[8:])
                    vcf["data"][key] = data_holder
            elif isPassUserFilters is not None:
                if isPassUserFilters(data_holder):
                    # passing_filters += 1
                    key = getKey(line, "chr_pos")
                    data_holder = parseFORMAT(data_holder, line[8:])
                    vcf["data"][key] = data_holder
                    vcf["data"][key]["vcf_line"] = line
            else:
                key = getKey(line, "chr_pos")
                data_holder = parseFORMAT(data_holder, line[8:])
                vcf["data"][key] = data_holder
                vcf["data"][key]["vcf_line"] = line
    
    # base = os.path.basename(path)
    # sample = os.path.splitext(base)[0]
    # logging.info("%s\t%s\t%s" % (sample, total_variant_count, passing_filters))
    
    # logging.debug(path)
    # for key in vcf["data"]:
    #     logging.debug(key)
    
    return vcf


