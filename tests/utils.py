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

from clinicalfilter.variant.snv import SNV
from clinicalfilter.variant.cnv import CNV
from clinicalfilter.trio_genotypes import TrioGenotypes

def make_var(chrom, pos,  alt, sex, info, keys, values, is_cnv=False, id='.',
        ref='A', filt='PASS', extra_info=None, format=None):
    ''' function to prevent duplication between create_snv and create_cnv
    '''
    
    if extra_info is not None:
        info = '{};{}'.format(info, extra_info)
    
    if format is not None:
        for (key, value) in format:
            keys += ':' + key
            values += ':' + value
    
    var = SNV
    if is_cnv:
        var = CNV
    
    return var(chrom, pos, id, ref, alt, filt, info=info, format=keys,
        sample=values, gender=sex)

def create_snv(sex, genotype, cq='missense_variant', hgnc='1001', chrom='1',
        pos='150', extra_info=None, format=None):
    ''' create a default SNV
    '''
    
    alt = 'G'
    
    info = 'HGNC_ID={0};CQ={1};DENOVO-SNP'.format(hgnc, cq)
    keys = 'GT:DP:TEAM29_FILTER:PP_DNM'
    values = '{}:50:PASS:0.99'.format(genotype)
    
    return make_var(chrom, pos, alt, sex, info, keys, values,
        extra_info=extra_info, format=format)

def create_cnv(sex, genotype, cq='missense_variant', hgnc='1001', chrom='1',
        pos='150', extra_info=None, format=None):
    ''' create a default CNV
    '''
    
    alt = '<DUP>'
    
    svlen = 5000
    info = 'CQ={0};HGNC_ID={1};HGNC_ALL={1};END={2};SVLEN={3}'.format(cq, hgnc, int(pos) + svlen, svlen)
    keys = 'INHERITANCE:DP'
    values = '{}:50'.format(genotype)
    
    return make_var(chrom, pos, alt, sex, info, keys, values, is_cnv=True,
        extra_info=extra_info, format=format)

def create_variant(child_sex, cq, hgnc, chrom='1', pos='150'):
    ''' create a default TrioGenotypes variant
    '''
    
    # generate a test variant
    child = create_snv(child_sex, '0/1', cq, hgnc, chrom)
    mom = create_snv('F', '0/0', cq, hgnc, chrom)
    dad = create_snv('M', '0/0', cq, hgnc, chrom)
    
    return TrioGenotypes(chrom, pos, child, mom, dad)

def make_vcf_header(sample_id='sample'):
    ''' generate a test VCF header
    '''
    
    lines = ['##fileformat=VCFv4.1\n',
        '##fileDate=2014-01-01\n',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(sample_id)]
    
    return lines

def make_vcf_line(chrom=1, pos=1, ref='G', alts='T',
        cq='missense_variant', genotype='0/1', extra=None):
    ''' generate a VCF line suitable for the unit tests
    
    Args:
        chrom: chromosome as string
        pos: nucleotide position of the variant
        ref: reference allele
        alts: comma-separated alternate alleles
        cq: vep consequence string. Can be '|' separated (for multiple
            genes) and/or ',' separated (for multiple alt alleles).
        genotype: SNV-based genotype for the variant e.g. '0/1'
        extra: additional fields for the INFO, e.g. 'HGNC=TEST'
    
    Returns:
        string for VCF line
    '''
    
    info = 'CQ={}'.format(cq)
    if extra is not None:
         info += ';' + extra
    
    return '{}\t{}\t.\t{}\t{}\t1000\tPASS\t{}\tGT:DP\t{}:50\n'.format(chrom,
        pos, ref, alts, info, genotype)

def make_minimal_vcf():
    ''' construct the bare minimum of lines for a VCF file
    '''
    
    variants = []
    variants.append(make_vcf_line(pos=100))
    variants.append(make_vcf_line(pos=200))
    
    return make_vcf_header() + variants
