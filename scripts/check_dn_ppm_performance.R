### Script to check the performance of changing the PP_DNM threshold for 
### putative de novos.

library(gdata)

# define some paths and folders
USER_DIR = "/nfs/users/nfs_j/jm33"
FREEZE_DIR = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/"

VALIDATED_FILE = file.path(FREEZE_DIR, "DNG_Validation_1133trios_20140130.tsv")
# VALIDATED_SHEET = "Exome variants reviewed"
REFERENCE_FILE = file.path(USER_DIR, "clinical_reporting.2014-05-01.txt")
VARIANTS_FILE = file.path(USER_DIR, "clinical_reporting.pp_dnm_0001.txt")


open_validated_variants <- function(filename){
    # open the list of reviewed variants, with their validated status
    
    # validated = read.xls(filename, sheet = sheet)
    validated = read.table(filename, sep = "\t", header = TRUE)
    
    return(validated)
}

open_variants <- function(filename) {
    # opens a dataframe of variants that pass clinical filtering
    
    tab = read.table(filename, sep = "\t", header = TRUE, 
        blank.lines.skip = TRUE)
    
    return(tab)
}

get_de_novos <- function(data) {
    # extract de novo variants from a dataframe
    
    # get the de novo SNVs from the autosomes and allosomes
    autosomal = data[data$chrom != "X" & data$trio_genotype == "1/0/0", ]
    male_allosomal = data[data$sex == "M" & data$chrom == "X" & 
        (data$trio_genotype == "2/0/0" | data$trio_genotype == "2/0/2"), ]
    female_allosomal = data[data$sex == "F" & data$chrom == "X" & 
        data$trio_genotype == "1/0/0", ]
    
    all = rbind(autosomal, male_allosomal, female_allosomal)
    
    return(all)
}

count_decisions <- function(de_novos, validated) {
    # count the de novos by their validated decisions
    
    # assign unique IDs, so we can get intersections and relative complements
    de_novos$unique_id = paste(de_novos$proband, de_novos$chrom, 
        de_novos$position, sep = "_")
    validated$unique_id = paste(validated$person_stable_id, validated$chr, 
        validated$pos, sep = "_")
    
    # count the number of de novos included in the filtered variants that have
    # been validated or not
    validated = validated[validated$unique_id %in% de_novos$unique_id, ]
    # print(table(de_novo_validated$DECISION))
    
    # find the de novos missing from the filtered variants, and count the number
    # of variants that have been validated or not
    passed = nrow(validated[validated$validation_result == "DNM", ])
    failed = nrow(validated[validated$validation_result != "DNM", ])
    # print(table(missing$DECISION))
    print(c(failed, passed / (passed + failed)))
     
}

main <- function() {
    # only open the validated variants file once (since it is slow to load)
    if (!exists("validated")) {
        validated <<- open_validated_variants(VALIDATED_FILE)
    }
    
    # reference = open_variants(REFERENCE_FILE)
    variants = open_variants(VARIANTS_FILE)
    
    # ref_de_novos = get_de_novos(reference)
    current_de_novos = get_de_novos(variants)
    
    count_decisions(current_de_novos, validated)
    
}

main()




