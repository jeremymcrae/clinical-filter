# code to check overlap or not between different sets of filtered variants

require(venneuler)
require(Cairo)
require(gdata)

HOME = "/nfs/users/nfs_j/jm33"
DIAGNOSES_PATH="/nfs/ddd0/Data/datafreeze/1133trios_20131218/Diagnosis_Summary_1133_20140328.xlsx"

exclude_cnv_based_variants <- function(initial) {
    # trim out CNVs, and variants which are present because they are compound hets
    # with a CNV
    cnvs = initial[grepl("DUP|DEL", initial$trio_genotype), c("proband", "chrom", "gene")]
    remove = cnvs[0, ]
    for (row_num in 1:nrow(cnvs)) {
        cnv = cnvs[row_num, ]
        cnv_genes = gsub(",", "|", cnv$gene)
        
        # get the variants for the same proband and on the same chrom
        matching = initial[initial$proband == cnv$proband & initial$chrom == cnv$chrom, ]
        
        # convert variant genes to vectors, and find which variants match any CNV gene
        genes = strsplit(matching$gene, ",")
        idx = sapply(genes, function(x) any(grepl(cnv_genes, x)))
        
        remove = rbind(remove, matching[idx, ])
    }
    
    # and drop the CNV related variants from the dataset
    initial = initial[!row.names(initial) %in% row.names(remove), ]
    
    return(initial)
}

# identifies the VCF path for a person from within a PED file
get_vcf_path <- function(sample_id, ped_path) {
    ped = read.table(ped_path, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    names(ped) = c("family_id", "individual_id", "father_id", "mother_id", "sex", "affected", "vcf")
    
    sample = ped[ped$individual_id == sample_id, ]
    
    return(sample$vcf)
}

# constructs a pedigree file for a single proband (and their parents)
construct_pedigree_file <- function(proband_id, original_path, new_path) {
    
    ped = read.table(original_path, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    names(ped) = c("family_id", "individual_id", "father_id", "mother_id", "sex", "affected", "vcf")
    
    proband_line = ped[ped$individual_id == proband_id, ]
    father_line = ped[ped$individual_id == proband_line$father_id, ]
    mother_line = ped[ped$individual_id == proband_line$mother_id, ]
    
    family = rbind(proband_line, father_line, mother_line)
    
    write.table(family, file=new_path, sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
}

# run clinical filtering on a set of variants, debugging each time, so that can
# identify the reason why each individual variant failed the clinical filtering
check_variant_fail <- function(variants, ped_path) {
    
    variants$fail_reason = NA
    for (row_num in 1:nrow(variants)) {
        proband_id = variants$proband[row_num]
        new_path = tempfile(fileext=".ped")
        construct_pedigree_file(proband_id, ped_path, new_path)
        
        pos = variants$position[row_num]
        chrom = variants$chrom[row_num]
        
        print(c(proband_id, chrom, pos))
        
        command = "python3" 
        args = c("/nfs/users/nfs_j/jm33/apps/clinical-filter/src/main/python/clinical_filter.py",
            "--ped",  new_path, 
            "--filter", "/nfs/users/nfs_j/jm33/apps/clinical-filter/config/filters.txt", 
            "--tags", "/nfs/users/nfs_j/jm33/apps/clinical-filter/config/tags.txt",
            "--deprecated-genes", "/nfs/users/nfs_j/jm33/apps/clinical-filter/config/ddg2p_deprecated_hgnc.txt", 
            "--output", "test.txt", 
            "--syndrome-regions", "/lustre/scratch113/projects/ddd/resources/decipher_syndrome_list_20140428.txt",
            "--known-genes", "/nfs/ddd0/Data/datafreeze/1133trios_20131218/DDG2P_with_genomic_coordinates_20131107_updated_TTN.tsv", 
            "--debug-chrom", chrom, 
            "--debug-pos", pos)
        
        reason = system2(command, args, stdout=TRUE)
        # allow for character(0) standard out
        if (is.character(reason) & length(reason) == 0L) { reason = NA}
        # trim multi-line output
        reason = reason[1]
        if (nchar(reason) > 100) { reason = substr(reason, 1, 100) }
        
        if (is.na(reason)) {
            path = get_vcf_path(proband_id, ped_path)
            string = paste("^", chrom, "\t", pos, sep="")
            grep_check = grep(string, readLines(path))
            if (is.integer(grep_check) & length(grep_check) == 0L) { reason = "variant not in VCF"}
        }
        print(reason)
        
        variants$fail_reason[row_num] = reason
        file.remove(new_path)
    }
    return(variants)
}

# aggregate the fail categories, so we can easily make tables
combine_fail_categories <- function(variants) {
    
    variants$fail_category = NA
    
    # I've manually set all the failure categories. Most of the variants fail for
    # the MAF, polyphen prediction, and VQSLOD categories.
    variants$fail_category[grepl("VQSLOD", variants$fail_reason)] = "failed VQSLOD"
    variants$fail_category[grepl("not smaller_than 0\\.01", variants$fail_reason)] = "failed MAF"
    variants$fail_category[grepl("low MAF in non-biallelic", variants$fail_reason)] = "failed MAF"
    variants$fail_category[grepl("variant not in VCF", variants$fail_reason)] = "variant not in VCF"
    variants$fail_category[grepl("polyphen prediction", variants$fail_reason)] = "polyphen changed"
    variants$fail_category[grepl("failed DENOVO-SNP\\/INDEL check", variants$fail_reason)] = "not called by denovogear"
    variants$fail_category[grepl("gatk_LowQual", variants$fail_reason)] = "failed gatk filter"
    variants$fail_category[grepl("gatk_StandardFilters", variants$fail_reason)] = "failed gatk filter"
    variants$fail_category[grepl("TEAM29_FILTER", variants$fail_reason)] = "failed denovogear filters"
    variants$fail_category[grepl("PP_DNM", variants$fail_reason)] = "failed PP_DNM"
    variants$fail_category[grepl("failed CQ", variants$fail_reason)] = "consequence changed"
    variants$fail_category[grepl("failed HGNC", variants$fail_reason)] = "gene symbol not picked up"
     variants$fail_category[is.na(variants$fail_reason)] = "likely missing compound pair"
    
    return(variants)
}

main <- function() {
    diagnoses = read.xls(DIAGNOSES_PATH, sheet="Exome variants reviewed", stringsAsFactors=FALSE)
    diagnoses = subset(diagnoses, select=c("proband", "chrom", "position", "DECISION"))
    initial_path = file.path(HOME, "clinical_reporting.2014-10-30.txt")
    current_path = file.path(HOME, "clinical_reporting.txt")
    
    initial = read.table(initial_path, sep="\t", stringsAsFactors=FALSE, header=TRUE, blank.lines.skip=TRUE)
    current = read.table(current_path, sep="\t", stringsAsFactors=FALSE, header=TRUE, blank.lines.skip=TRUE)
    initial = exclude_cnv_based_variants(initial)
    
    # standardise the individuals in the datasets
    current = current[current$proband %in% initial$proband, ]
    initial = initial[initial$proband %in% current$proband, ]
    
    both = merge(initial, current, by=c("proband", "chrom", "position", "sex"), all=TRUE)
    both = merge(both, diagnoses, by=c("proband", "chrom", "position"), all.x=TRUE)
    
    initial_only = both[is.na(both$ref.alt_alleles.y), ]
    current_only = both[is.na(both$ref.alt_alleles.x), ]
    
    # show the proportions for each trio genotype type
    table(initial$trio_genotype)/nrow(initial)
    table(current$trio_genotype)/nrow(current)
    
    # check the reasons why each variant failed
    initial_only = check_variant_fail(initial_only, "~/samples.ped")
    current_only = check_variant_fail(current_only, "~/exome_reporting.ped")
    
    # categorize the fail reasons
    initial_only = combine_fail_categories(initial_only)
    current_only = combine_fail_categories(current_only)
    
    # get the numbers of variants in each fail category (in the inital_only
    # cross-tabulate by reporting decision)
    table(initial_only$fail_category, initial_only$DECISION)
    table(current_only$fail_category)
    
    # plot a venn diagram of the numbers of variants in each group
    v = venneuler(c(previous=nrow(initial_only), current=nrow(current_only), "previous&current"=(nrow(both)-(nrow(initial_only) + nrow(current_only)))))
    Cairo(file="variants_venn_diagram.pdf", type="pdf", height=10, width=10, units="cm")
    plot(v)
    dev.off()
    
    # construct_pedigree_file("DDDP100149", "~/samples.ped", "temp.ped")
}


  


