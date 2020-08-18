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

from __future__ import unicode_literals

import io
import json
import re

def get_header_positions(file_handle, columns):
    """ get a dictionary of column positions from a header line
    
    Args:
        file_handle: file handle for known genes file.
        columns: list of column names eg ["gene", "start", "stop"].
    
    Returns:
        dictionary of column positions indexed by the column name.
    """
    
    line = file_handle.readline().strip().split("\t")
    
    positions = {}
    for column in columns:
        positions[column] = line.index(column)
    
    return positions

def parse_gene_line(line, header):
    """ adds a gene to the known genes dictionary by parsing the data line
    
    Args:
        line: list of columns for gene line
        header: dictionary of column positions indexed by column name
    
    Returns:
        HGNC ID and dictionary entry for gene
    """
    
    hgnc_id = line[header["hgnc_id"]]
    symbol = line[header["gene"]]
    status = line[header["type"]]
    inheritance = line[header["mode"]]
    mechanism = line[header["mech"]]
    
    gene = {}
    if re.search(r',', inheritance):
        inh_split = inheritance.split(',')
        gene['inh'] = {}
        for inher in inh_split:
            gene['inh'][inher.title()] = set([mechanism])
    else:
        gene['inh'] = {inheritance: set([mechanism])}
    gene["symbol"] = symbol
    gene["status"] = set([status.lower()])
    gene["start"] = int(line[header["start"]])
    gene["end"] = int(line[header["stop"]])
    gene["chrom"] = line[header["chr"]]
    
    # some genes are listed with an inheritance mode of "Both", which means
    # the gene has been observed in disorders with both monoallelic and
    # biallelic inheritance. Make sure the monoallelic and biallelic modes
    # are shown for the gene.
    if inheritance == "Both":
        gene["inh"]["Monoallelic"] = set([mechanism])
        gene["inh"]["Biallelic"] = set([mechanism])

    return hgnc_id, gene

def open_known_genes(path):
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
    """
    
    if path is None:
        return None
    
    # only include genes with sufficient DDG2P status
    allowed = set(["confirmed dd gene", "probable dd gene", "both rd and if"])
    
    known = {}
    with io.open(path, "r", encoding="latin_1") as handle:
        # get the positions of the columns in the list of header labels
        columns = ["gene", "type", "mode", "mech", "start", "stop", "chr", "hgnc_id"]
        header = get_header_positions(handle, columns)
        
        for line in handle:
            line = line.strip().split("\t")
            hgnc_id, gene = parse_gene_line(line, header)
            
            if len(gene['status'] & allowed) == 0:
                continue
            
            if hgnc_id not in known:
                known[hgnc_id] = gene
            else:
                for mode in gene['inh']:
                    if mode not in known[hgnc_id]['inh']:
                        known[hgnc_id]['inh'][mode] = set()
                    known[hgnc_id]['inh'][mode] |= gene['inh'][mode]
                
                # merge entries for genes with multiple modes or mechanisms
                known[hgnc_id]['status'] |= gene['status']
    
    if len(known) == 0:
        raise ValueError("No genes found in the file, check the line endings")

    return known

def open_cnv_regions(path):
    """ opens a file listing CNV regions
    
    Args:
        path: path to CNV regions file
    
    Returns:
        dictionary of copy number values, indexed by (chrom, start end) tuples
    """
    
    if path is None:
        return None
    
    cnv_regions = {}
    with open(path) as handle:
        header = handle.readline()
        for line in handle:
            line = line.strip().split("\t")
            
            chrom = line[5]
            start = line[3]
            end = line[4]
            copy_number = line[2]
            
            key = (chrom, start, end)
            cnv_regions[key] = copy_number
    
    return cnv_regions

def open_last_base_sites(path):
    ''' open a set of sites at the last base of an exon which can potentially
    alter the consequence to a LoF consequence.
    
    Args:
        path: path to last base sites file, or None
    
    Returns:
        Set of sites as (chrom, pos) tuples. Can be empty set if path is None.
    '''
    
    if path is None:
        return set([])
    
    with open(path) as handle:
        last_base = json.load(handle)
        # convert everything to tuples, which are hashable into a set
        last_base = set([ (x[0], int(x[1])) for x in last_base ])
    
    return last_base

def open_x_lr2_file(path):
    '''open file containing sum of mean log 2 ratios for each proband.
    Args:
        path: path to x_lr2 file
    Returns:
        Set of proband and sum xl2r as a dict
     '''

    if path is None:
        return {}

    sumxlr2 = {}
    with open(path) as handle:
        for line in handle:
            linedata = line.split()
            if linedata[0].startswith('DDD'):
                sumxlr2[linedata[0]] = linedata[1]

    return sumxlr2
            
