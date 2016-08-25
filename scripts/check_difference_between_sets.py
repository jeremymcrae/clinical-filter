# code to check overlap or not between different sets of filtered variants

import os
import subprocess
import tempfile

import pandas

HOME = "/nfs/users/nfs_j/jm33/clinical_reporting_output"

def get_vcf_path(sample_id, ped_path):
    ''' identifies the VCF path for a person from within a PED file
    '''
    ped = read.table(ped_path)
    sample = ped[ped['individual_id'] == sample_id].squeeze()
    
    return sample['vcf']

def construct_pedigree_file(proband_id, original_path, new_path):
    ''' constructs a pedigree file for a single proband (and their parents)
    '''
    
    ped = read.table(original_path)
    
    proband_line = ped[ped['individual_id'] == proband_id]
    father_line = ped[ped['individual_id'] == proband_line['father_id']]
    mother_line = ped[ped['individual_id'] == proband_line['mother_id']]
    
    family = rbind(proband_line, father_line, mother_line)
    
    family.to_csv(new_path, sep="\t", index=False)

def check_variant_fail(variant, ped_path):
    ''' run clinical filtering on a set of variants, debugging each time, so that can
    identify the reason why each individual variant failed the clinical filtering
    
    Args:
     variant: data frame or list for single variant, providing proband,
         chrom and position info
     ped_path: path to a pedigree file that contains family information for
         the proband and the probands parents. Note that the ped file can contain
         info for other families, since we construct a ped file specifically for
         the proband of interest within this function, using the info from the
         ped_path argument.
    
    Examples:
        check_variant_fail(list(proband="DDDP100099", chrom=6, position=152712420),
            "~/samples.ped")
    '''
    
    proband_id = variant["proband"]
    temp_ped = tempfile.NamedTemporaryFile(suffix=".ped")
    construct_pedigree_file(proband_id, ped_path, temp_ped.name)
    
    pos = variant["position"]
    chrom = variant["chrom"]
    
    command = ['python3', "/nfs/users/nfs_j/jm33/apps/clinical-filter/bin/clinical_filter.py",
        "--ped",  temp_ped.name,
        "--output", "test.txt",
        "--syndrome-regions", "/lustre/scratch113/projects/ddd/resources/decipher_syndrome_list_20140428.txt",
        "--known-genes", "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2014-11-04/DDG2P/dd_genes_for_clinical_filter",
        "--debug-chrom", chrom,
        "--debug-pos", pos]
    
    reason = subprocess.check_output(command)
    # allow for character(0) standard out
    # if (is.character(reason) & length(reason) == 0L) { reason = NA}
    # # trim multi-line output
    reason = reason.split('\n')[0]
    if len(reason) > 100:
        reason = reason[0:100]
    
    # if reason == '':
    #     path = get_vcf_path(proband_id, ped_path)
    #     string = paste("^", chrom, "\t", pos, sep="")
    #     grep_check = grep(string, readLines(path))
    #     if (is.integer(grep_check) & length(grep_check) == 0L) { reason = "variant not in VCF"}
    
    print(proband_id, chrom, pos, reason)
    temp_ped.close()
    
    return reason

def combine_fail_categories(variants):
    ''' aggregate the fail categories, so we can easily make tables
    '''
    
    category = pandas.Series([None] * len(variants))
    
    # I've manually set all the failure categories. Most of the variants fail for
    # the MAF, polyphen prediction, and VQSLOD categories.
    category[variants['reason'].contains("LOW_VQSLOD")] = "failed LOW_VQSLOD"
    category[variants['reason'].contains("not smaller_than 0\\.01")] = "failed MAF"
    category[variants['reason'].contains("low MAF in non-biallelic")] = "failed MAF"
    category[variants['reason'].contains("variant not in VCF")] = "variant not in VCF"
    category[variants['reason'].contains("polyphen prediction")] = "polyphen changed"
    category[variants['reason'].contains("failed DENOVO-SNP\\/INDEL check")] = "not called by denovogear"
    category[variants['reason'].contains("gatk_LowQual")] = "failed gatk filter"
    category[variants['reason'].contains("gatk_StandardFilters")] = "failed gatk filter"
    category[variants['reason'].contains("TEAM29_FILTER")] = "failed denovogear filters"
    category[variants['reason'].contains("PP_DNM")] = "failed PP_DNM"
    category[variants['reason'].contains("failed CQ")] = "consequence changed"
    category[variants['reason'].contains("failed HGNC")] = "gene symbol not picked up"
    category[variants['reason'].isnull()] = "likely missing compound pair"
    
    return category

def get_difference(table_1, table_2):
    '''
    '''
    
    table1 = table_1.copy()
    table2 = table_2.copy()
    
    # code the positions as integers (some floats can get in otherwise)
    table1['position'] = table1['position'].astype(int)
    table2['position'] = table2['position'].astype(int)
    
    table1['key'] = table1["proband"] + '-' + table1["chrom"] + '-' + \
        table1['position'].map(str) + '-' + table1['sex']
    
    table2['key'] = table2["proband"] + '-' + table2["chrom"] + '-' + \
        table2['position'].map(str) + '-' + table2['sex']
    
    return table_1[~table1['key'].isin(table2['key'])]

def main():
    initial_path = os.path.join(HOME, "clinical_reporting.2016-08-24.unmodified.txt")
    current_path = os.path.join(HOME, "clinical_reporting.2016-08-25.overdominance.txt")
    
    initial = pandas.read_table(initial_path)
    current = pandas.read_table(current_path)
    
    initial = initial[~initial['proband'].isnull()]
    current = current[~current['proband'].isnull()]
    
    initial_only = get_difference(initial, current)
    current_only = get_difference(current, initial)
    
    # show the proportions for each trio genotype type
    initial['trio_genotype'].value_counts()/len(initial)
    current['trio_genotype'].value_counts()/len(current)
    
    basename = os.path.basename(current_path)
    initial_only.to_csv(basename + '.removed.txt', sep='\t', index=False)
    current_only.to_csv(basename + '.added.txt', sep='\t', index=False)
    
    # # check the reasons why each variant failed
    # initial_only$reason = apply(initial_only, 1, check_variant_fail, CURRENT_PED)
    # current_only$reason = apply(current_only, 1, check_variant_fail, ORIGINAL_PED)
    #
    # # categorize the fail reasons
    # initial_only$category = combine_fail_categories(initial_only)
    # current_only$category = combine_fail_categories(current_only)

if __name__ == '__main__':
    main()
