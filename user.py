'''
Read user filter and grouping files and return them as python dict

10 Dec 2011
sa9@sanger.ac.uk
'''

import os
import sys
import logging

def parseGenes(path='DDGP-reportable.txt'):
    """Loads list of known disease causative genes.
    
    We obtain a liost of genes that are known to be involved in disorders, so that we can screen 
    variants for being in these known genes, which makes them better candidates for being causative.
    We currently only use genes that are monoallelic or biallelic.
    
    Args:
        path: path to tab-separated file listing known disease-causing genes.
    
    Returns:
        A dictionary of genes, so we can check variants for inclusion in the set. The dictionary is 
        indexed by gene ID to the corresponding inheritance value
    
    Raises:
        IOError: an error when the gene file path is not specified correctly.
    """
    
    if not os.path.exists(path):
        raise IOError("Path to known gene file does not exist")
    
    known_genes = {}
    f = open(path, 'r')
    
    # allow for gene files with different column names and positions
    header = f.readline().strip().split("\t")
    if "DDG2P_Status" in header:
        gene_label = "Gene"
        confirmed_status_label = "DDG2P_Status"
        inheritance_label = "Inheritance"
    elif "type" in header:
        gene_label = "gene"
        confirmed_status_label = "type"
        inheritance_label = "mode"
    else:
        raise ValueError("The gene file doesn't contain any expected header column names")
    
    # get the positions of the columns in the list of header labels
    gene_column = header.index(gene_label)
    confirmed_status_column = header.index(confirmed_status_label)
    inheritance_column = header.index(inheritance_label)
    
    # only include genes with sufficient DDG2P status
    allowed_confirmed_statuses = ["Confirmed DD Gene", "Probable DD gene", "Both DD and IF"]
    
    for line in f:
        line = line.strip().split("\t")
        gene_ID = line[gene_column]
        gene_confirmed_status = line[confirmed_status_column]
        gene_inheritance = line[inheritance_column]
        
        # ignore genes with insufficient evidence
        if gene_confirmed_status not in allowed_confirmed_statuses:
            continue 
        
        if gene_ID not in known_genes:
            known_genes[gene_ID] = set([])
        
        known_genes[gene_ID].add(gene_inheritance)
        
        # some genes are listed with an inheritance mode of "Both", which means the gene has been
        # observed in disorders with both monoallelic and biallelic inheritance. Make sure the
        # monoallelic and biallelic modes are shown for the gene.
        if gene_inheritance == "Both":
            known_genes[gene_ID].add("Monoallelic")
            known_genes[gene_ID].add("Biallelic")
        
        # # currently only report mono and biallelic genes
        # if gene_inheritance == 'Biallelic':
        #     known_genes[gene_ID] = 2
        # elif gene_inheritance == "Monoallelic":
        #     known_genes[gene_ID] = 1
    
    if len(known_genes) == 0:
        raise ValueError("No genes found in the file, check the line endings")
    
    return known_genes

def create_person_ID_mapper(path="/nfs/ddd0/Data/datafreeze/1139trios_20131030/person_sanger_decipher.private.txt"):
    """creates a dictionary of IDs to map between different ID systems.
    
    We occasionally need to convert between different ID schemas (eg between DDD person IDs and 
    DECIPHER person IDs). 
    
    Args:
        path: path to tab-separated file listing the alternate IDs
    
    Returns:
        dictionary with current ID and alternate ID as key value pairs for different individuals.
    """
    if not os.path.exists(path):
        raise IOError("Path to known ID converter file does not exist")
    
    ID_converter = {}
    f = open(path, "r")
    for line in f:
        line = line.strip().split("\t")
        person_ID = line[0]
        alternate_ID = line[1]
        
        ID_converter[person_ID] = alternate_ID
    
    return ID_converter

def parseFilters(path):
    """Opens user-defined criteria for filtering variants.

    Open a tab-separated text file with VCF criteria for variants (such as only including function 
    altering nonsynonymous variants, or only including variants where any of the 1000 Genomes 
    cohorts have minor allele frequencies less than 0.01). Each line in the file is for a separate
    criteria, with the columns being filter_label (corresponding to a VCF INFO category), 
    filter_type (eg list, greater_than, smaller_than), and the actual criteria (could be a comma-
    separated list, or a numeric value).
    
    Args:
        path: path to text file defining the filters.

    Returns: 
        A dictionary of VCF INFO categories with corresponding criteria. For example: 
        
        {'VCQ': '(list', ['ESSENTIAL_SPLICE_SITE', 'STOP_GAINED']), 
        'AF_MAX': ('smaller_than', 0.01)}

    Raises:
        IOError: An error when the filter file path is not specified correctly.
        ValueError: An error when apparently numeric values cannot be converted to floats.
    """

    # check that the file exists
    if not os.path.exists(path):
        raise IOError("path to filter file does not exist")

    # open the path and load the filter criteria
    mydict = {}
    f = open(path, 'r')
    for line in f:
        if not line.startswith("#"):
            label, condition, values = line.strip().split("\t")
            
            # split the comma-separated strings into sets
            if condition == "list":
                values = set(values.split(","))
            
            # split nested strings into nested lists
            if condition == "multiple_not":
                values = values.strip("()").split("), (")
                for position in range(len(values)):
                    values[position] = tuple(values[position].split(", "))
            
            # convert numeric values to floats
            if condition in set(["greater_than", "smaller_than", "equal"]):
                try:
                    values = float(values)
                except ValueError:
                    print "Please check your filter file or correct this value. Here is the full record:"
                    print line
                    raise ValueError("One of these values couldn't be converted to float. Filter (%s) VCF value (%s)" % (condition, values))

            # convert numeric lists to lists of floats
            if condition == 'range':
                try:
                    values = [float(x) for x in values]
                except ValueError:
                    print "Please check your filter file or correct this value. Here is the full record:"
                    print line
                    raise ValueError("One of these values couldn't be converted to float. Filter (%s) VCF value (%s)" % (condition, values))

            mydict[label] = (condition, values)
    f.close()
    return mydict

def parseWeights(path, order=False):
    """Opens user-defined heirarchy of weights for prioritising variants.

    Opens a tab-separated text file with user-defined weights, which are used to prioritise 
    candidate variants for further examination. 

    Args:
        path: path to the text file containing the weights.
        order: whether to return the weights as a list from parseOrders().

    Returns:
        A dictionary of VFC INFO categories with corresponding criteria. For example: 

        {'AF_MAX': {'smaller_than': [(0.01, 1, 'RARE_1KG')]}, 
        'VCQ': {'list': [('ESSENTIAL_SPLICE_SITE,STOP_GAINED', 1, 'LOF&FUNC'), 
        ('SYNONYMOUS_CODING', 0, 'SILN')]}}

    Raises:
        IOError: An error when the weights file is not specified correctly.

    NOTE: This function appears similar to parseFilters(). Why is it necessary? Is it because there 
    are more criteria listed her, or that we get to specify an order for the weight list?
    """

    # check that the weights file exists
    if not os.path.exists(path):
        raise IOError("path to weights file does not exist")

    if order == False:
        mydict = {}
        f = open(path,'r')
        for line in f:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                tag = line[0]
                condition = line[1]
                weight = float(line[3])
                label = line[4]

                values = line[2]
                # sometimes the values are a comma-separated list, split them into a list
                if values.count(",") > 0:
                    values = values.split(",")
                
                # convert numeric values to floats
                if condition in set(["greater_than", "smaller_than", "equal"]):
                    try:
                        values = float(values)
                    except ValueError:
                        print "Please check your filter file or correct this value. Here is the full record:"
                        print line
                        raise ValueError("One of these values couldn't be converted to float. Filter (%s) VCF value (%s)" % (condition, values))

                # convert numeric lists to lists of floats
                if condition == 'range':
                    try:
                        values = [float(x) for x in values]
                    except ValueError:
                        print "Please check your filter file or correct this value. Here is the full record:"
                        print line
                        raise ValueError("One of these values couldn't be converted to float. Filter (%s) VCF value (%s)" % (condition, values))
               
                # multiple conditions can apply to the same tag, account for those as sub-dictionaries
                if not mydict.has_key(tag):
                    mydict[tag] = {condition:[(values, weight, label)]}
                else:
                    if not mydict[tag].has_key(condition):
                        mydict[tag][condition] = [(values, weight, label)]
                    else:
                        mydict[tag][condition].append((values, weight, label))
        f.close()
        #mydict = addEmptyClasses(mydict)
        return mydict
    else:
        # the parseOrders function gets a list of the tags used in the weights file.
        # This might be simpler as:
        # if order == True:
        #     return list(mydict.keys())
        # or even passing this list from earlier in EVAR_trios.py For instance, changing:
        # self.orders  = user.parseOrders(columns_path, weights_path))
        # to:
        # self.orders  = user.parseOrders(columns_path, list(self.weights.keys()))
        weightsOrders = []
        f = open(path,'r')
        for line in f:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                if line[0] not in weightsOrders:
                    weightsOrders.append(line[0])
        f.close()
        return weightsOrders


def parseOrders(orderPath, weightsPath):
    """Gets a list of which VCF attributes to include in the output.

    Opens a tab-separated text file (typically named output.txt) and obtains a dictionary of the 
    attributes listed in there. The primary bit seems to be a line starting with 'output', which 
    has a comma-separated list of attributes. Also includes the weights file attributes (which 
    should be included in the output file, so I don't understand why this is needed yet).

    Args:
        orderPath: path to the text file containing the output attributes
        weightsPath: path to the weights path, for additional attributes

    Returns:
        A dictionary of VCF attributes included in the output and weights: For example:

        {'output': ['CHROM', 'POS', 'ID', 'REF'...], 
        'weights': ['AF_MAX', 'VCQ']}

    Raises:
        IOError: An error when the orderPath is not specified correctly.
    """

    if not os.path.exists(orderPath):
        raise IOError("path to output file does not exist")

    mydict = {}
    f = open(orderPath,'r')
    for line in f:
        if not line.startswith("#"):
            line = line.strip().split("\t")
            key = line[0]
            
            # split the comma-separated values into a list if they exist
            try:
                values = line[1].split(",")
            except:
                values = line[1]
            mydict[key] = values
    f.close()
    
    # it may be easier to pass the weights dict keys to this function, rather than generating from
    # the weights file
    mydict["weights"] = parseWeights(weightsPath, True)
    
    return mydict

def parseTags(path):
    """Opens alternate identifiers for values in nonstandard VCF files.

    The EVA_report needs to know what the GT and GN tags are called in the VCF file. Depending on 
    where the VCF files has originated from, these can be named differently. This function opens a 
    tab-separated file, and gets values to map one nomenclature to another. The GT is expected in 
    the 9 column of the VCF file and its value in the 10th column. The gene ID (GN) is expected in 
    INFO column (7th column).

    Args:
        path: path to tags file (typically named tags.txt)

    Returns:
        Strings for each of the desired tags eg CQ_tag = "CQ", or GN_tag = "VGN".
    
    Raises:
        IOError: An error when the tags path is not specified correctly.
    """

    if not os.path.exists(path):
        raise IOError("path to tags file does not exist")

    GN_tag = "gene"
    CQ_tag = "consequence"
    MAF_tag = "MAX_MAF"
    GT_tag = "genotype"
    
    tags_dict = {GN_tag: "", CQ_tag: "", MAF_tag: "", GT_tag: ""}
    f = open(path,'r')
    for line in f:
        if line.startswith("#"):
            continue
        
        line = line.strip().split("\t")
        key = line[0]
        value = line[1]
        
        # sometimes the identifier matches a list of possible tags
        if "," in value:
            value = value.split(",")
        
        # if there is only a single value, we might as well insert it into a list, for consistency
        if type(value) == str:
            value = [value]
        
        tags_dict[key] = value
        
    return tags_dict

def printFileContent(path):
    """ prints the text content of a file.
    """
    f = open(path, 'r')
    for line in f:
        if not line.startswith('#'):
            logging.info("#" + line.strip())
    f.close()

def getFileContent(path):
    """Returns the text content of a file.
    """
    f = open(path, 'r')
    text = f.read()
    f.close()
    return text

def checkUserInput(path, type):
    '''
    Check if filters.txt , hierarchy.txt , output.txt and tags.txt 
    have the right format otherwise exit and warn the user.
    '''
    pass
    #TODO
