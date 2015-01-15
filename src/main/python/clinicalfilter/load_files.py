""" loads configuration files for clinical-filter.py
"""

import os
import io
import time
import sys

def open_file(path):
    """ opens a file handle, allowing for permission errors
    
    Occasionally the filesystem can block opening a file, if another process is
    accessing it. We simply wait a few seconds before trying to open the file
    again (but cut out after a few retries, in case something else has gone
    wrong).
    
    Args:
        path: path to a file
    
    Returns:
        a file object handle, or raises an error
    """
    
    retry = 0
    
    while retry < 5:
        try:
            return io.open(path, "r", encoding="latin_1")
        except PermissionError as e:
            # if we get an error, wait a few seconds for other processes to
            # release the file, before retrying
            time.sleep(3)
            retry += 1
    
    raise IOError("cannot access file, even after " + str(retry) + " tries.\n" \
        + "This is often seen when multiple processes try to access a file,\n" \
        + "and the lustre filesystem is being stressed.")

def open_known_genes(path="DDGP-reportable.txt"):
    """Loads list of known disease causative genes.
    
    We obtain a list of genes that are known to be involved in disorders, so
    that we can screen variants for being in these known genes, which makes them
    better candidates for being causative.
    
    Args:
        path: path to tab-separated file listing known disease-causing genes.
    
    Returns:
        A dictionary of genes, so we can check variants for inclusion in the 
        set. The dictionary is indexed by gene ID to the corresponding 
        inheritance value.
    
    Raises:
        IOError: an error when the gene file path is not specified correctly.
    """
    
    f = open_file(path)
    
    # allow for gene files with different column names and positions
    header = f.readline().strip().split("\t")
    
    # get the positions of the columns in the list of header labels
    gene_column = header.index("gene")
    confirmed_status_column = header.index("type")
    inheritance_column = header.index("mode")
    mechanism_column = header.index("mech")
    start_column = header.index("start")
    stop_column = header.index("stop")
    chrom_column = header.index("chr")
    
    # only include genes with sufficient DDG2P status
    allowed_statuses = set(["Confirmed DD Gene", "Probable DD gene", \
        "Both DD and IF"])
    
    excluded_genes = set([])
    known_genes = {}
    for line in f:
        line = line.strip().split("\t")
        gene = line[gene_column]
        confirmed_status = line[confirmed_status_column]
        inheritance = line[inheritance_column]
        mechanism = line[mechanism_column]
        
        # ignore genes with insufficient evidence
        if confirmed_status not in allowed_statuses:
            excluded_genes.add(gene)
            continue 
        
        if gene not in known_genes:
            known_genes[gene] = {"inh": {}, "status": set()}
        
        if inheritance not in known_genes[gene]["inh"]:
            known_genes[gene]["inh"][inheritance] = set()
        
        known_genes[gene]["inh"][inheritance].add(mechanism)
        known_genes[gene]["status"].add(confirmed_status)
        known_genes[gene]["start"] = int(line[start_column])
        known_genes[gene]["end"] = int(line[stop_column])
        known_genes[gene]["chrom"] = line[chrom_column]
        
        # some genes are listed with an inheritance mode of "Both", which means 
        # the gene has been observed in disorders with both monoallelic and 
        # biallelic inheritance. Make sure the monoallelic and biallelic modes 
        # are shown for the gene.
        if inheritance == "Both":
            if "Monoallelic" not in known_genes[gene]["inh"]:
                known_genes[gene]["inh"]["Monoallelic"] = set()
            if "Biallelic" not in known_genes[gene]["inh"]:
                known_genes[gene]["inh"]["Biallelic"] = set()
            known_genes[gene]["inh"]["Monoallelic"].add(mechanism)
            known_genes[gene]["inh"]["Biallelic"].add(mechanism)
    
    if len(known_genes) == 0:
        raise ValueError("No genes found in the file, check the line endings")
    
    f.close()
    
    return known_genes, excluded_genes

def create_person_ID_mapper(path="/nfs/ddd0/Data/datafreeze/1139trios_20131030/person_sanger_decipher.private.txt"):
    """creates a dictionary of IDs to map between different ID systems.
    
    We occasionally need to convert between different ID schemas (eg between
    DDD person IDs and DECIPHER person IDs). 
    
    Args:
        path: path to tab-separated file listing the alternate IDs
    
    Returns:
        dictionary with current ID and alternate ID as key value pairs for 
        different individuals.
    """
    
    converter = {}
    with open(path) as f:
        for line in f:
            # don't bother to include the maternal or paternal IDs, since we 
            # only want the probands IDs
            if ":mat" in line or ":pat" in line:
                continue
            
            line = line.split("\t")
            person_id = line[0]
            alternate_id = line[1]
            
            if person_id not in converter:
                converter[person_id] = alternate_id
    
    return converter

def open_filters(path):
    """ Opens user-defined criteria for filtering variants.

    Open a tab-separated text file with VCF criteria for variants (such as only 
    including function altering nonsynonymous variants, or only including 
    variants where any of the 1000 Genomes cohorts have minor allele frequencies
    less than 0.01). Each line in the file is for a separate criteria, with the
    columns being filter_label (corresponding to a VCF INFO category), 
    filter_type (eg list, greater_than, smaller_than), and the actual criteria 
    (could be a comma-separated list, or a numeric value).
    
    Args:
        path: path to text file defining the filters.

    Returns: 
        A dictionary of VCF INFO categories with corresponding criteria, e.g.
        
        {"VCQ": "(list", ["ESSENTIAL_SPLICE_SITE", "STOP_GAINED"]), 
        "AF_MAX": ("smaller_than", 0.01)}

    Raises:
        IOError: An error when the filter file path is not specified correctly.
        ValueError: An error when apparently numeric values cannot be converted
            to floats.
    """
    
    f = open_file(path)

    # open the path and load the filter criteria
    filters = {}
    for line in f:
        if line.startswith("#"):
            continue
            
        label, condition, values = line.strip().split("\t")
        
        # split the comma-separated strings into sets
        if condition == "list":
            values = set(values.split(","))
        
        # convert numeric values to floats
        elif condition in set(["greater_than", "smaller_than", "equal"]):
            values = float(values)
        
        # convert numeric lists to lists of floats
        elif condition == "range":
            values = [float(x) for x in values]

        filters[label] = (condition, values)
            
    f.close()
    
    return filters

def open_tags(path):
    """ Opens alternate identifiers for values in nonstandard VCF files.

    The report needs to know what the GT and GN tags are called in the VCF file.
    Depending on where the VCF files has originated from, these can be named 
    differently. This function opens a tab-separated file, and gets values to 
    map one nomenclature to another. The GT is expected in the 9 column of the 
    VCF file and its value in the 10th column. The gene ID (GN) is expected in 
    INFO column (7th column).

    Args:
        path: path to tags file (typically named tags.txt)

    Returns:
        dictionary listing possible values for each of the desired tags eg 
        CQ_tag = "CQ", or GN_tag = "VGN"
    """
    
    f = open_file(path)
    
    GN_tag = "gene"
    CQ_tag = "consequence"
    MAF_tag = "MAX_MAF"
    GT_tag = "genotype"
    
    tags_dict = {GN_tag: "", CQ_tag: "", MAF_tag: "", GT_tag: ""}
    for line in f:
        if line.startswith("#"):
            continue
        
        line = line.strip().split("\t")
        key = line[0]
        value = line[1]
        
        # sometimes the identifier matches a list of possible tags
        if "," in value:
            value = value.split(",")
        
        # if there is only a single value, we might as well insert it into a 
        # list, for consistency
        if type(value) == str:
            value = [value]
        
        tags_dict[key] = value
        
    return tags_dict

def open_cnv_regions(path):
    """ opens a file listing CNV regions
    
    Args:
        path: path to CNV regions file
    
    Returns:
        dictionary of copy number values, indexed by (chrom, start end) tuples
    """
    
    f = open_file(path)
    header = f.readline().strip().split("\t")
    
    cnv_regions = {}
    for line in f:
        line = line.strip().split("\t")
        
        chrom = line[5]
        start = line[3]
        end = line[4]
        copy_number = line[2]
        
        key = (chrom, start, end)
        cnv_regions[key] = copy_number
    
    f.close()
    
    return cnv_regions
