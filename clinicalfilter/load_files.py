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

def add_gene_to_known_genes(line, header, known_genes):
    """ adds a gene to the known genes dictionary by parsing the data line
    
    Args:
        line: list of columns for gene line
        header: dictionary of column positions indexed by column name
        known_genes: dictionary of known genes
    
    Returns:
        known genes, with new gene added
    """
    
    gene = line[header["gene"]]
    status = line[header["type"]]
    inheritance = line[header["mode"]]
    mechanism = line[header["mech"]]
    
    known_genes[gene]["inh"][inheritance].add(mechanism)
    known_genes[gene]["status"].add(status)
    known_genes[gene]["start"] = int(line[header["start"]])
    known_genes[gene]["end"] = int(line[header["stop"]])
    known_genes[gene]["chrom"] = line[header["chr"]]
    
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
    
    return known_genes

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
    """
    
    # only include genes with sufficient DDG2P status
    allowed = set(["Confirmed DD Gene", "Probable DD gene", "Both DD and IF"])
    
    known_genes = {}
    with io.open(path, "r", encoding="latin_1") as handle:
        # get the positions of the columns in the list of header labels
        columns = ["gene", "type", "mode", "mech", "start", "stop", "chr"]
        header = get_header_positions(handle, columns)
        
        for line in handle:
            line = line.strip().split("\t")
            gene = line[header["gene"]]
            status = line[header["type"]]
            inheritance = line[header["mode"]]
            mechanism = line[header["mech"]]
            
            # ignore genes with insufficient evidence
            if status not in allowed:
                continue
            
            if gene not in known_genes:
                known_genes[gene] = {"inh": {}, "status": set()}
            
            if inheritance not in known_genes[gene]["inh"]:
                known_genes[gene]["inh"][inheritance] = set()
            
            known_genes = add_gene_to_known_genes(line, header, known_genes)
    
    if len(known_genes) == 0:
        raise ValueError("No genes found in the file, check the line endings")
    
    return known_genes

def create_person_ID_mapper(path):
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
    with open(path) as handle:
        for line in handle:
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

def open_cnv_regions(path):
    """ opens a file listing CNV regions
    
    Args:
        path: path to CNV regions file
    
    Returns:
        dictionary of copy number values, indexed by (chrom, start end) tuples
    """
    
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
