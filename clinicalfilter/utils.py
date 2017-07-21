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


import os
import io
import sys
import gzip
import hashlib

from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV

IS_PYTHON3 = sys.version_info.major == 3

def open_vcf(path):
    """ Gets a file object for an individual's VCF file.
    
    Args:
        path: path to VCF file (gzipped or text format).
        
    Returns:
        A file handle for the VCF file.
    """
    
    if not os.path.exists(path):
        raise OSError("VCF file not found at: " + path)
    
    extension = os.path.splitext(path)[1]
    
    if extension == ".gz":
        # python2 gzip opens in text, but same mode in python3 opens as
        # bytes, avoid with platform specific code
        handle = gzip.open(path, "r")
        if IS_PYTHON3:
            handle = gzip.open(path, "rt")
    elif extension in [".vcf", ".txt"]:
        handle = io.open(path, "r", encoding="latin_1")
    else:
        raise OSError("unsupported filetype: " + path)
    
    return handle

def get_vcf_header(path):
    """ Get the header lines from a VCF file.
    
    Args:
        path: path to VCF file, or file handle.
    
    Returns:
        a list of lines that start with "#", which are the header lines.
    """
    
    is_handle = False
    try:
        vcf = open_vcf(path)
    except TypeError:
        vcf = path
        is_handle = True
    
    current_pos = vcf.tell()
    vcf.seek(0)
    
    header = []
    for line in vcf:
        if not line.startswith("#"):
            break
        
        header.append(line)
    
    vcf.seek(current_pos)
    
    # this is a bit awkward, but if we've passed in a path, we want to close
    # the is_handle, otherwise we leave an opened file in unit tests.
    if not is_handle:
        vcf.close()
    
    return header

def exclude_header(vcf):
    """ removes the header from a VCF file object
    
    We remove the header from the VCF file, since the header is ~200 lines
    long, and an exome VCF file is 100,000 lines long, so it's better to
    remove the header once, rather than continually check if lines are part
    of the header as we traverse the VCF. We simply run through the VCF
    until we find a non-header line, then seek back to the start of that
    line.
    
    Args:
        vcf: handler for a VCF file
    """
    
    current_pos = vcf.tell()
    
    while vcf.readline().startswith("#"):
        current_pos = vcf.tell()
    
    vcf.seek(current_pos)

def construct_variant(line, gender, mnvs=None):
    """ constructs a Variant object for a VCF line, specific to the variant type
    
    Args:
        line: list of elements of a single sample VCF line:
            [chrom, position, snp_id, ref_allele, alt_allele, quality,
            filter_value, info, format_keys, format_values]
        gender: gender of the individual to whom the variant line belongs
            (eg "1" or "M" for male, "2", or "F" for female).
    
    Returns:
        returns a Variant object
    """
    
    chrom, pos, var_id, ref, alt, qual, filter_val, info, format_keys, sample = line[:10]
    
    # CNVs are found by their alt_allele values, as either <DUP>, or <DEL>
    Var = SNV
    if alt in ["<DUP>", "<DEL>"]:
        Var = CNV
    
    mnv_code = None
    if mnvs is not None and (chrom, int(pos)) in mnvs:
        mnv_code = mnvs[(chrom, int(pos))]
    
    var = Var(chrom, pos, var_id, ref, alt, qual, filter_val, info, format_keys,
        sample, gender, mnv_code=mnv_code)
    
    if var.is_cnv():
        var.fix_gene_IDs()
    
    return var

def get_vcf_provenance(person):
    """ get provenance information for a VCF
    
    Args:
        person: Person object for an individual, or None if the person doesn't exist
    
    Returns:
        returns a tuple of sha1 VCF file hash, name of VCF file (without
        directory), and date the VCF file was generated
    """
    
    if person is None:
        return ('NA', 'NA', 'NA')
    
    path = person.get_path()
    
    # get the SHA1 hash of the VCF file (in a memory efficient manner)
    BLOCKSIZE = 65536
    checksum = hashlib.sha1()
    with open(path, "rb") as handle:
        buf = handle.read(BLOCKSIZE)
        while len(buf) > 0:
            checksum.update(buf)
            buf = handle.read(BLOCKSIZE)
    
    checksum = checksum.hexdigest()
    basename = os.path.basename(path)
    
    date = None
    for line in get_vcf_header(path):
        if line.startswith("##fileDate"):
            date = line.strip().split("=")[1]
            break
    
    # some VCF files lack the fileDate in the header, get it from the path
    if date is None:
        date = os.path.splitext(basename)[0]
        date = date.split(".")[2]
    
    return (checksum, basename, date)
