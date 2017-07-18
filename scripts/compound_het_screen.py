'''
Copyright (c) 2017 Genome Research Ltd.

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

import argparse
import tempfile
import gzip
from copy import deepcopy
from itertools import chain

from clinicalfilter.utils import open_vcf, get_vcf_header, exclude_header, \
    construct_variant

import pysam

def get_options():
    ''' parse command line options
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--id', help='ID of individual (must be in the VCF)')
    parser.add_argument('--vcf', help='path to VCF file')
    parser.add_argument('--bam', help='path to BAM file')
    parser.add_argument('--out', help='path to write filtered VCF to')
    
    return parser.parse_args()

def remove_phased_hets(person, vcf_path, bam_path, output_vcf_path):
    ''' screen out putative compound hets in phase from a VCF.
    
    We identify putative compound hets when we lack parental data, but sequence
    reads can identify candidates where both sites are in the same read, which
    means both sites are inherited from the same parent, and excludes the site
    as being a compound het.
    
    Args:
        person: sample ID for person, which
        vcf_path: path to probands VCF
        bam_path: path to probands BAM sequence
        output_vcf_path: path to write filtered VCF to
    '''
    
    phased = [ x for x in get_compound(vcf_path, person) if in_phase(bam_path, x) ]
    phased = set([ x for sublist in phased for x in sublist ])
    
    initial_vcf = open_vcf(vcf_path)
    header = get_vcf_header(initial_vcf)
    exclude_header(initial_vcf)
    
    output_vcf = gzip.open(output_vcf_path, 'wt')
    output_vcf.writelines(header)
    
    for line in initial_vcf:
        record = construct_variant(line.split('\t'), 'F')
        # only write out variants which are not 'compound hets' in phase
        key = (record.chrom, record.position, record.ref_allele, record.alt_alleles)
        if key not in phased:
            output_vcf.write(line)

def get_compound(vcf_path, sample_id):
    ''' pull out the compound hets, grouped by gene
    
    Args:
        vcf_path: path to VCF
        sample_id: sample ID for individual in VCF.
    
    Returns:
        list of lists of (chrom, pos, ref, alts) tuples for the variants in a
        compound het.
    '''
    
    vcf = open_vcf(vcf_path)
    exclude_header(vcf)
    
    genes = {}
    for line in vcf:
        line = line.split('\t')
        variant = construct_variant(line, 'F')
        variant.add_format(line[8], line[9])
        
        if 'compound_het' not in variant.info['ClinicalFilterType']:
            continue
        
        # only check sites in singletons, which always have inheritance=unknown
        if variant.format['INHERITANCE'] != 'unknown':
            continue
        
        for symbol in variant.info['ClinicalFilterReportableHGNC']:
            if symbol not in genes:
                genes[symbol] = []
            
            genes[symbol].append((variant.chrom, variant.position,
                variant.ref_allele, variant.alt_alleles))
    
    return genes.values()

def parse_cigar(read):
    ''' get the positions of indels within the read
    
    We give the indels positions that correspond to the reference genome
    coordinates.
    '''
    
    codes = {0: 'BAM_CMATCH', 1: 'BAM_CINS', 2: 'BAM_CDEL', 3: 'BAM_CREF_SKIP',
        4: 'BAM_CSOFT_CLIP', 5: 'BAM_CHARD_CLIP', 6: 'BAM_CPAD',
        7: 'BAM_CEQUAL', 8: 'BAM_CDIFF'}
    
    indels = {}
    pos = read.reference_start
    for (key, length) in read.cigar:
        key = codes[key]
        if key == 'BAM_CINS':
            indels[pos - 1] = length
        elif key == 'BAM_CDEL':
            indels[pos] = -length
        
        if key != 'BAM_CINS':
            pos += length
    
    return indels

def read_has_alt(read, pos, ref, alt):
    ''' checks whether the sequence at a site is for a 'ref' or 'alt' allele
    
    Args:
        read: AlignmentRead from pileup
        ref: sequence for reference allele
        alt: sequence for alternate alllele
    
    Returns:
        'ref' or 'alt' depending on whether the allele is for the reference
        sequence, or alternate to the reference sequence.
    '''
    
    if alt == '*':
        alt = ''
    
    is_indel = len(ref) != 1 or len(alt) != 1
    
    if is_indel:
        length = (len(alt) - len(ref))
        indels = parse_cigar(read)
        
        if pos in indels:
            return length == indels[pos]
        return False
    
    try:
        offset = read.get_reference_positions().index(pos - 1)
    except ValueError:
        # if we try to check a SNV for a position which has a deleyion covering
        # the site of the SNV, then the SNV position won't exists within the
        # reference positions. This won't match the alt allele, since we check
        # for indels elsewhere.
        return False
    
    # if the variant is for a SNV, get the base call as a string
    seq = read.query_alignment_sequence[offset:offset + 1]
    return seq == alt

def check(reads, variant):
    for read in reads:
        if read is None:
            return False
        
        if read.reference_start < variant[1] < read.reference_end:
            return read_has_alt(read, variant[1], variant[2], variant[3][0])
    
    return None

def end(pos, ref, alt):
    ''' get the end position of a variant
    
    Args:
        pos: integer for nucleotide position
        ref: reference allele
    
    Returns:
        end position of variant as integer
    '''
    
    delta = abs(len(ref) - len(alt))
    
    return pos + delta - 1

def get_mate(bam, read):
    ''' get the mate pair of a read (or None if no mate)
    
    Args:
        bam: pysam.AlignmentFile for bam sequence
        read: pysam.AlignedSegment for sequence read
    
    Returns:
        AlignedSegment sequence read, or None if no paired read available
    '''
    
    try:
        return bam.mate(read)
    except ValueError:
        return None

def check_overlap(pos, ref, alt, others):
    ''' check if a variant has any overlapping alternates
    
    Args:
        pos: nucleotide position on chromosome
        ref: reference allele
        alt: alternate allele
        others: list of (pos, ref, alt) for the compound het partners
    '''
    # account for variants where two variants overlap. These can have
    # different alt representations to the VCF, so won't be picked up by
    # looking for matching alts in the read.
    return all([ (end(x[1], x[2], x[3]) >= pos >= x[1]) or \
        (end(x[1], x[2], x[3]) >= end(pos, ref, alt) >= x[1]) for x in others ])

def in_phase(bam_path, variants):
    ''' determine if the phase data suggests the gene has a true compound het
    
    for all of the variants in the list, check that there is at least one
    partner variant where either the phase differs or is unknown
    
    Args:
        bam: path to bam file for proband
        variants: list of (chrom, pos, ref, alts) tuples for variants in a
            compound het.
    
    Returns:
         boolean for whether the variants phases data suggests the compound hets
         are in phase.
    '''
    
    bam = pysam.AlignmentFile(bam_path)
    
    phasing = []
    for (chrom, pos, ref, alts) in sorted(variants):
        
        others = deepcopy(variants)
        others.remove((chrom, pos, ref, alts))
        
        alt = alts[0]
        
        pileup = bam.pileup(chrom, pos-1, pos, truncate=True, stepper='nofilter')
        reads = ( read.alignment for x in pileup for read in x.pileups )
        mates = ( get_mate(bam, x) for x in reads )
        
        counts = {'in_phase': 0, 'not_in_phase': 0}
        for (read, mate) in zip(reads, mates):
            
            if check_overlap(pos, ref, alt, others):
                counts['in_phase'] += 1
                continue
            
            if not read_has_alt(read, pos, ref, alt):
                continue
            
            possible = [ check([read, mate], x) for x in others ]
            if False in possible or None in possible:
                counts['not_in_phase'] += 1
            else:
                counts['in_phase'] += 1
        
        phasing.append(counts)
    
    # Many variants lack phase data, due to being too far apart for the reads to
    # provide evidence. We'll only exclude variants where we can confidently
    # classify the variants as being in phase.
    
    # exclude variants with all zero counts, due to distant variants
    if all(( x['in_phase'] == 0 and x['not_in_phase'] == 0 for x in phasing )):
        return False
    
    return sum([ x['in_phase'] == 0 for x in phasing ]) < 2

def get_header(vcf):
    ''' extract the header from a VCF, so we can use it to write a new VCF
    
    The reason this is a function rather than accessing the header attribute
    is some VCFs have a broken header, where they lack an INFO field description
    for CANDIDATE_MNV. Otherwise the code segfaults when writing variants with
    CANDIDATE_MNV fields. I have amended upstream code to define CANDIDATE_MNV
    in the header in the output VCF, but I want to cover VCFs already produced.
    
    Returns:
        header from VCF.
    '''
    
    header = vcf.header
    
    mnv_header = '##INFO=<ID=CANDIDATE_MNV,Number=.,Type=String,' \
        'Description="Code for candidate multinucleotide variants. Field is ' \
        'only included if the translated MNV differs from both of the SNV ' \
        'translations. There are five possibilities: alternate_residue_mnv=' \
        'MNV translates to a residue not in SNVs, masked_stop_gain_mnv=' \
        'MNV masks a stop gain, modified_stop_gained_mnv=MNV introduces a ' \
        'stop gain, modified_synonymous_mnv=MNV reverts to synonymous, ' \
        'modified_protein_altering_mnv=synonymous SNVs but missense MNV."\n'
    
    if 'CANDIDATE_MNV' not in header.info.keys():
        header.add_line(mnv_header)
    
    return header

def main():
    args = get_options()
    remove_phased_hets(args.id, args.vcf, args.bam, args.out)

if __name__ == '__main__':
    main()
