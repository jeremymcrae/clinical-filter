'''
Copyright (c) 2016 Wellcome Trust Sanger Institute

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

import re
from collections import namedtuple

import tabix

from clinicalfilter.load_vcfs import LoadVCFs

coding_cq = set(["transcript_ablation", "splice_donor_variant",
    "splice_acceptor_variant", "stop_gained", "frameshift_variant",
    "stop_lost", "start_lost", "initiator_codon_variant",
    "inframe_insertion", "inframe_deletion", "missense_variant",
    "protein_altering_variant", "transcript_amplification",
    "splice_region_variant", "incomplete_terminal_codon_variant",
    "synonymous_variant", "stop_retained_variant",
    "coding_sequence_variant"])
    
AA_CODE = {"AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
    "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
    "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F"}

def get_mnv_candidates(path):
    ''' identify MNV candidates, and their MNV consequences within a VCF.
    
    Args:
        path: path to VCF
    
    Returns:
        list of (variant, mnv_consequence) tuples, where variant is (chrom, pos)
    '''
    
    loader = LoadVCFs(total_trios=0, known_genes=None, last_base=None,
        debug_chrom=None, debug_pos=None)
    with loader.open_vcf_file(path) as vcf:
        loader.exclude_header(vcf)
        header = loader.get_vcf_header(vcf)
        pairs = find_nearby_variants(vcf)
    
    # ensure variants are not indels, are coding, and pairs alter the same amino acid
    vcf = tabix.open(path)
    pairs = screen_pairs(vcf, pairs, is_not_indel)
    pairs = screen_pairs(vcf, pairs, is_coding)
    pairs = same_aa(vcf, pairs)
    
    pattern = re.compile('[ACGT]')
    
    candidates = []
    for pair in pairs:
        var1, var2 = list(get_matches(vcf, pair))
        cq = check_mnv_consequence(var1, var2, pattern)
        candidates.append((pair[0], cq))
        candidates.append((pair[1], cq))
    
    return candidates

def find_nearby_variants(vcf, threshold=2):
    ''' find variants in close proximity, regardless of allele or consequence
    
    Args:
        vcf: pysam.VariantFile object, for the current VCF.
        threshold: distance in base-pairs for variants to be nearby.
    
    Returns:
        list of chromosome pairs (each member as a (chrom, pos) tuple)
    '''
    
    # This simply looks for VCF lines adjacent to each other where the variants
    # in close proximity. I tried checking in a window of 100 preceeding lines,
    # but this didn't net any extra MNV candidates, so just checking the
    # previous line should be fine.
    nearby = []
    previous = ('0', -10000)
    for variant in vcf:
        
        # chrom = variant.chrom
        # pos = variant.pos
        chrom, pos, rest = variant.split('\t', 2)
        pos = int(pos)
        
        if chrom != previous[0]:
            previous = (chrom, pos)
            continue
        
        delta = abs(previous[1] - pos)
        # make sure the match isn't for the same variant, which avoids a bug
        # with duplicate lines
        if delta < threshold and previous[1] != pos:
            nearby.append([previous, (chrom, pos)])
        
        previous = (chrom, pos)
    
    return nearby

def parse_vcf_line(line, Variant):
    ''' parse a VCF line into a useable form. This loosly mimics the pysam setup
    
    Args:
        line: variant line from VCF, already split by tabs.
        Variant: namedtuple for a variant. Omits the samples and format, since
            we don't need them for identifying MNVs.
    
    Returns:
        namedtuple, with the correct entries filled in from the VCF line
    '''
    
    chrom, pos, id, ref, alts, qual, filter, info_str = line[0:8]
    pos = int(pos)
    alts = alts.split(',')
    
    # create a dictionary of values from the info fields
    info = {}
    for item in info_str.split(";"):
        try:
            key, value = item.split("=")
        except ValueError:
            key, value = item, True
        info[key] = value
    
    return Variant(chrom, pos, id, ref, alts, qual, filter, info)

def get_matches(vcf, pair):
    ''' find VCF lines matching a pair of coordinate tuples
    
    Args:
        vcf: pytabix file for VCF
        pair: a list of (chrom, pos) tuples
    
    Yields:
        VariantRecord for matching variants
    '''
    
    Variant = namedtuple('Variant', ['chrom', 'pos', 'id', 'ref', 'alts', 'qual',
        'filter', 'info'])
    
    # define the coordinates
    chrom = pair[0][0]
    positions = set([ x[1] for x in pair ])
    start = min(positions)
    end = max(positions)
    
    # pull out the matching VCF variant entries
    # for var in vcf.fetch(chrom, start-1, end):
    for var in vcf.query(chrom, start-1, end):
        var = parse_vcf_line(var, Variant)
        if var.pos in positions:
            yield var

def is_not_indel(var):
    ''' check if a variant is not an indel
    
    Args:
        var: VariantRecord for a single variant
    '''
    ref_snv = len(var.ref) == 1
    alt_snv = any([ len(y) == 1 and y != '*' for y in var.alts ])
    
    return ref_snv and alt_snv

def is_coding(var):
    ''' check if any of the consequence annotations are coding
    
    Args:
        var: VariantRecord for a single variant
    '''
    cq = set()
    for x in var.info['CQ'].split(','):
        cq |= set(x.split('|'))
    
    return len(cq & coding_cq) > 0

def screen_pairs(vcf, pairs, func):
    ''' exclude proximal pairs where at least one partner fails.
    
    This is a generic function, and we pass in a function to check each variant
    
    Args:
        vcf: pysam.cbcf.VariantFile for VCF
        pairs: list of chromosome pairs (each member as a (chrom, pos) tuple)
        func: screening function e.g. is_coding or is_not_indel
    
    Returns:
        subset of the list of pairs, where the pairs have to pass.
    '''
    
    cleaned = []
    for pair in pairs:
        checks = [ func(x) for x in get_matches(vcf, pair) ]
        
        assert len(checks) == 2
        
        if all(checks):
            cleaned.append(pair)
    
    return cleaned

def same_aa(vcf, pairs):
    ''' exclude proximal pairs where the partners are not in the same amino acid
    '''
    
    cleaned = []
    for pair in pairs:
        
        aa = []
        for x in get_matches(vcf, pair):
            # splice_region variants can be outside CDS, make these fail
            if 'Protein_position' not in x.info:
                aa += [1, 2]
            else:
                aa.append(x.info['Protein_position'])
        
        if len(set(aa)) == 1:
            cleaned.append(pair)
    
    return cleaned

def translate(seq):
    ''' translate a codon sequence
    '''
    
    seq = seq.upper()
    
    protein = ''
    n = 3
    for codon in (seq[i:i+n] for i in range(0, len(seq), n)):
        protein += AA_CODE[codon]
    
    return protein

def get_codons(var1, var2, pattern):
    ''' get the refernce, and variant codons for the SNV/MNV codons
    
    Args:
        var1: VariantRecord for a single variant
        var2: VariantRecord for a single variant
        pattern: compiled regex pattern for uppercases bases
    
    Returns:
        dictionary of sequences for the reference, both SNVs and the combined
        MNV codon.
    '''
    
    alt_1 = var1.alts[0]
    alt_2 = var2.alts[0]
    
    # the codon sequences are in the info 'Codons' field, which can include
    # multiple transcripts ('|' separared). So far, alternate transcripts were
    # missing codon sequences (coded as '.'). Trim these out by removing the '.'
    # value from the codon sequences set. This won't allow different transcripts
    # with differing codon positions, but that shouldn't happen.
    codons = set(var1.info['Codons'].split(',')[0].split('|')) - set(['.'])
    codons2 = set(var2.info['Codons'].split(',')[0].split('|')) - set(['.'])
    
    # check we only have a single codon/transcript entry per variant
    assert len(codons) == 1
    assert len(codons2) == 1
    
    codons = list(codons)[0].split('/')
    codons2 = list(codons2)[0].split('/')
    
    # get the individual variant codons
    reference = codons[0]
    var1 = codons[1]
    var2 = codons2[1]
    
    # get the positions of the modified bases, which is shown by the capitalised
    # base in the codon sequence e.g. 'tGt' or 'Agg'.
    var1_pos = pattern.search(var1, 0).start()
    var2_pos = pattern.search(var2, 0).start()
    
    # figure out the position of the unmodified codon base
    unmodified_pos = set([0, 1, 2]) - set([var1_pos, var2_pos])
    unmodified_pos = list(unmodified_pos)[0]
    
    bases = dict(zip([unmodified_pos, var1_pos, var2_pos],
        [reference, var1, var2]))
    
    # create a merged codon that includes both variant bases
    mnv = ''
    for key in sorted(bases):
        mnv += bases[key][key]
    
    return {'reference': reference, 'snv1': var1, 'snv2': var2, 'mnv': mnv}

def check_mnv_consequence(var1, var2, pattern):
    ''' get the multinucleotide variant consequence vs the individual SNVs.
    
    we've got several scenarios:
    1. all the codons produce the same amino acid - ie there is no difference
       between any of the codons.
    2. the MNV codon gives the same change as one or both of of the SNV codons.
    3. the MNV codon gives a different change compared to the SNV codons.
       3a. the MNV codon alters the protein, where the SNV codons do not.
       3b. the MNV codon is synonymous, where the SNV codons are not.
       3c. the MNV codon is stop gained, where the SNV codons are not.
       3d. the MNV prevents a stop_gained from one or both of the SNV codons.
       3e. the MNV codon and the SNV codons produce different residues.
    
    Args:
        var1: VariantRecord for a single variant
        var2: VariantRecord for a single variant
        pattern: compiled regex pattern for uppercases bases
    
    Returns:
        string indicating the MNV change conpared to the SNV codons
    '''
    
    codons = get_codons(var1, var2, pattern)
    ref = translate(codons['reference'])
    snv1 = translate(codons['snv1'])
    snv2 = translate(codons['snv2'])
    mnv = translate(codons['mnv'])
    
    if set([ref, snv1, snv2, mnv]) == set([ref]):
        change = 'unmodified_synonymous_mnv'
    elif mnv in [snv1, snv2]:
        change = 'unmodified_protein_altering_mnv'
    else:
        if ref == snv1 == snv2 and ref != mnv:
            change = 'modified_protein_altering_mnv'
        elif ref == mnv and (ref != snv1 or start_aa != snv2):
            change = "modified_synonymous_mnv"
        elif mnv == '*':
            change = 'modified_stop_gained_mnv'
        elif '*' in [snv1, snv2]:
            change = 'masked_stop_gain_mnv'
        else:
            change = 'alternate_residue_mnv'
    
    return change
