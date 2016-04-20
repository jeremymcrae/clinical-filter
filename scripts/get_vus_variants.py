""" prioritise VUS variants from DDG2P/non-DDG2P variant lists.

Uses variants filtered by clinical filtering script. As some background, 
variants in DDG2P genes are prioritised for clinical review. The clinical 
filtering script can also output variants from all genes (not just DDG2P genes).
The caveat with the all-genes output is that we don't know the likely mode of 
inheritance. We prioritise variants from the all-gene analysis by only 
including:
    - validated functional de novo variants
    - rare loss-of-function homozygous, compound heterozygous and hemizygous 
      variants. 

"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import argparse
import os

LOF = set(["transcript_ablation","splice_donor_variant", \
    "splice_acceptor_variant", "frameshift_variant", "stop_gained", \
    "coding_sequence_variant"])

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="extract variants for DDD research track.")
    parser.add_argument("--vus", dest="vus_variants", help="VUS variant filename")
    parser.add_argument("-o", "--out", dest="output", help="output filename")
    
    args = parser.parse_args()
    
    return args.vus_variants, args.output

def open_variants(filename):
    """ open variants from a file, indexed by (ID, chrom, pos) tuples
    """
    
    f = open(filename, "r")
    
    variants = {}
    for line in f:
        # ignore header and blank lines
        if line == "\n" or line.startswith("proband"):
            continue
        
        # get the info for a unique key to identify the variant
        line = line.strip().split("\t")
        individual_id = line[0]
        sex = line[2]
        chrom = line[3]
        position = line[4]
        
        key = (individual_id, chrom, position)
        if individual_id not in variants:
            variants[individual_id] = {}
            variants[individual_id]['sex'] = sex
        
        variants[individual_id][key] = line
    
    return variants

def is_compound_het_lof(key, variants, lof_consequences):
    """ checks that both variants in a compound het are loss-of-function
    """
    
    line = variants[key]
    consequence = line[8].split(",")[0]
    inheritance = line[11]
    result = line[15]
    gene = line[5]
    
    sample_id = key[0]
    chrom = key[1]
    
    lof_compound_hets = 0
    for alt_key in variants:
        alt_sample_id = alt_key[0]
        alt_chrom = alt_key[1]
        
        line = variants[alt_key]
        alt_gene = line[5]
        
        alt_consequence = line[8].split(",")[0]
        alt_inheritance = line[11]
        alt_result = line[15]
        
        # skip past variant for different individuals, or on different chroms
        if alt_sample_id != sample_id or alt_chrom != chrom or alt_gene != gene:
            continue
        
        if alt_consequence in lof_consequences and "compound_het" in alt_result:
            lof_compound_hets += 1
    
    return lof_compound_hets >= 2

def get_lof_recessive_variants(variants, lof_consequences):
    """ gets recessive variants with loss-of-function consequences
    
    Args:
        variants: dict of variant lines, indexed by (ID, chrom, pos) tuples
        
    Returns:
        dictionary of LOF variants, indexed by (ID, chrom, pos) tuples
    """
    
    dominant_modes = set(["Monoallelic", "X-linked dominant"])
    
    lof_variants = []
    for person_id in variants:
        for key in variants[person_id]:
            line = variants[person_id][key]
            consequence = line[8].split(",")[0]
            inheritance = line[11]
            result = line[15]
            
            chrom = key[1]
            
            # if the variant is part of a compound het, both variants have to be
            # lof, so drop compound het variants where one of the pair is non-lof
            if "compound_het" in result and not is_compound_het_lof(key, variants[person_id], lof_consequences):
                if result == "compound_het":
                    continue
                line[15] = "single_variant"
            
            if consequence in lof_consequences and inheritance not in dominant_modes:
                lof_variants.append('\t'.join(line) + '\n')
    
    return lof_variants

def get_maternal_chrx_lof_in_males(variants, lof_consequences):
    '''
    '''
    
    mode = 'Hemizygous'
    
    lof_variants = []
    for person_id in variants:
        if variants[person_id]['sex'] not in ['M', 'Male', 'male']:
            continue
        
        for key in variants[person_id]:
            chrom = key[1]
            if chrom not in ['X', '23', 'chrX', 'chr23']:
                continue
            
            line = variants[person_id][key]
            consequence = line[8].split(",")[0]
            inheritance = line[11]
            
            if consequence in lof_consequences and mode in inheritance.split(','):
                lof_variants.append('\t'.join(line) + '\n')
    
    return lof_variants

def write_output(filename, variants):
    """ writes the output to a file
    """
    
    f = open(filename, "w")
    
    header = "proband\talternate_ID\tsex\tchrom\tposition\tgene\tmutation_ID"+ \
        "\ttranscript\tconsequence\tref/alt_alleles\tMAX_MAF\tinheritance" + \
        "\ttrio_genotype\tmom_aff\tdad_aff\tresult\tomim_match\n"
    f.write(header)
    
    f.write(variants)
    
    f.close()

def main():
    """
    """
    
    vus_file, output_file = get_options()
    
    # open the datasets for the script
    vus_vars = open_variants(vus_file)
    
    recessive = get_lof_recessive_variants(vus_vars, LOF)
    maternal_chrx = get_maternal_chrx_lof_in_males(vus_vars, LOF)
    
    variants = recessive + maternal_chrx
    
    write_output(output_file, variants)


if __name__ == '__main__':
    main()



    
