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

This script also:
    - annotates whether the gene is in the OMIM morbid map genes 
        (see http://compbio.charite.de/hudson/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/genes_to_diseases.txt)
    - filters out variants in genes with multiple (>1) common (MAF >= 0.05) LOF
      variants (as determined from 1000 Genomes datasets).
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import

import argparse
import os


def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description="examine mutation clustering \
        in genes")
    parser.add_argument("--vus", dest="vus_variants", help="VUS variant filename")
    parser.add_argument("--ddg2p", dest="ddg2p_variants", help="DDG2P variant \
        filename")
    parser.add_argument("--omim", dest="omim_morbid", help="OMIM morbid map \
        gene filename")
    parser.add_argument("--common-lofs", dest="common_lof", help="File listing \
        genes containing common (MAF > 0.05) LOF variants")
    parser.add_argument("--de-novos", dest="de_novos", help="File listing \
        validated de novo variants in probands")
    parser.add_argument("-o", "--out", dest="output", help="output filename")
    
    args = parser.parse_args()
    
    return args.vus_variants, args.ddg2p_variants, args.de_novos, \
        args.omim_morbid, args.common_lof, args.output

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
        chrom = line[3]
        position = line[4]
        
        key = (individual_id, chrom, position)
        if individual_id not in variants:
            variants[individual_id] = {}
        
        variants[individual_id][key] = line
    
    return variants

def open_omim_genes(filename):
    """ open a list of OMIM morbid map genes
    """
    
    f = open(filename, "r")
    
    omim = set([])
    for line in f:
        gene = line.strip()
        omim.add(gene)
    
    return omim

def open_common_lof(filename):
    """ open a list of genes containing common loss-of-function variants
    """
    
    f = open(filename, "r")
    
    # ignore header
    f.readline()
    
    genes = set([])
    for line in f:
        line = line.strip().split("\t")
        gene = line[0]
        count = line[1]
        
        if int(count) > 1:
            genes.add(gene)
    
    return genes

def open_validated_denovos(filename):
    """ gets (ID, chrom, pos) tuple keys for validated denovos
    """
    
    f = open(filename, "rU")
    
    validated = set([])
    for line in f:
        line = line.strip().split("\t")
        
        proband_id = line[0].split("_")[0]
        chrom = line[3]
        pos = line[4]
        validation_status = line[8]
        
        if validation_status == "DNM":
            validated.add((proband_id, chrom, pos))
    
    return validated

def remove_common_lof_genes(variants, lof_genes):
    """ trims variants in genes with common loss-of-function variants
    """
    
    without_common_lofs = {}
    for person_id in variants:
        for key in variants[person_id]:
            line = variants[person_id][key]
            hgnc_symbols = line[5].split(",")
            
            # print(hgnc_symbols)
            
            # find which genes for a variant match the common LOF genes
            common_lof_match = False
            for hgnc in hgnc_symbols:
                if hgnc in lof_genes:
                    common_lof_match = True
                    break
            
            if not common_lof_match:
                if person_id not in without_common_lofs:
                    without_common_lofs[person_id] = {}
                without_common_lofs[person_id][key] = line
    
    return without_common_lofs

def remove_ddg_overlap(all_vars, ddg2p_vars):
    """ removes all genes variants that overlap with the DDG2P variants
    """
    
    unique = {}
    for person_id in all_vars:
        if person_id not in ddg2p_vars:
            unique[person_id] = all_vars[person_id]
        else:
            for key in all_vars[person_id]:
                if key not in ddg2p_vars[person_id]:
                    if person_id not in unique:
                        unique[person_id] = {}
                    
                    unique[person_id][key] = all_vars[person_id][key]
    
    return unique

def get_functional_de_novos(variants, de_novo_keys):
    """ finds variants which are validated de novos
    """
    
    de_novos = {}
    for person_id in variants:
        for key in variants[person_id]:
            if key in de_novo_keys:
                if person_id not in de_novos:
                    de_novos[person_id] = {}
                de_novos[person_id][key] = variants[person_id][key]
        
    return de_novos

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
        
    print(key, lof_compound_hets >= 2)
    
    return lof_compound_hets >= 2

def get_lof_recessive_variants(variants):
    """ gets recessive variants with loss-of-function consequences
    
    Args:
        variants: dict of variant lines, indexed by (ID, chrom, pos) tuples
        
    Returns:
        dictionary of LOF variants, indexed by (ID, chrom, pos) tuples
    """
    
    lof_consequences = set(["transcript_ablation","splice_donor_variant", \
        "splice_acceptor_variant", "frameshift_variant", "stop_gained", \
        "coding_sequence_variant"])
    dominant_inheritances = set(["Monoallelic", "X-linked dominant"])
    
    lof_variants = {}
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
                print(key)
                if result == "compound_het":
                    continue
                line[15] = "single_variant"
            
            # print(key, consequence in lof_consequences, inheritance not in dominant_inheritances)
            
            if consequence in lof_consequences and inheritance not in dominant_inheritances:
                if person_id not in lof_variants:
                    lof_variants[person_id] = {}
                
                lof_variants[person_id][key] = line
    
    return lof_variants

def annotate_omim_gene(variants, omim_genes):
    """ annotate variants with the OMIM morbid map gene matches
    """
    
    annotated_vars = {}
    for person_id in variants:
        for key in variants[person_id]:
            line = variants[person_id][key]
            hgnc_symbols = line[5].split(",")
            
            # find which genes for a variant match the OMIM morbid map genes
            omim_match = []
            for hgnc in hgnc_symbols:
                if hgnc in omim_genes:
                    omim_match.append(hgnc)
            
            # annotate the variant with the omim match (or lack)
            if len(omim_match) == 0:
                omim_match = ["NA"]
            else:
                omim_match = [",".join(omim_match)]
            
            line = "\t".join(line + omim_match) + "\n"
            annotated_vars[key] = line
    
    return annotated_vars  

def write_output(filename, variants):
    """ writes the output to a file
    """
    
    f = open(filename, "w")
    
    header = "proband\talternate_ID\tsex\tchrom\tposition\tgene\tmutation_ID"+ \
        "\ttranscript\tconsequence\tref/alt_alleles\tMAX_MAF\tinheritance" + \
        "\ttrio_genotype\tmom_aff\tdad_aff\tresult\tomim_match\n"
    f.write(header)
    
    for key in sorted(variants):
        f.write(variants[key])
    
    f.close()

def main():
    """
    """
    
    vus_file, ddg2p_file, de_novo_file, omim_file, common_lof_file, output_file = get_options()
    
    # open the datasets for the script
    vus_vars = open_variants(vus_file)
    ddg2p_vars = open_variants(ddg2p_file)
    common_lof_genes = open_common_lof(common_lof_file)
    
    de_novos = open_validated_denovos(de_novo_file)
    de_novo_variants = get_functional_de_novos(vus_vars, de_novos)
    
    # filter the VUS variants
    vus_vars = remove_ddg_overlap(vus_vars, ddg2p_vars)
    vus_vars = get_lof_recessive_variants(vus_vars)
    vus_vars = remove_common_lof_genes(vus_vars, common_lof_genes)
    
    # put all the validated de novos back in
    for person_id in de_novo_variants:
        if person_id not in vus_vars:
            vus_vars[person_id] = {}
        vus_vars[person_id].update(de_novo_variants[person_id])
    
    # annotate the variants with whether their genes are OMIM morbid map genes, 
    # and write the variants to a file
    omim_genes = open_omim_genes(omim_file)
    vus_vars = annotate_omim_gene(vus_vars, omim_genes)
    # print(vus_vars.keys())
    write_output(output_file, vus_vars)


if __name__ == '__main__':
    main()



    