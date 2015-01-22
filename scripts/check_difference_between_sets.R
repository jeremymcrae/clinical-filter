# code to check overlap or not between different sets of filtered variants

require(venneuler)
require(Cairo)
require(gdata)

HOME = "/nfs/users/nfs_j/jm33"
DIAGNOSES_PATH="/nfs/ddd0/Data/datafreeze/1133trios_20131218/Diagnosis_Summary_1133_20140328.xlsx"
ORIGINAL_PED = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/family_relationships.private.txt"
CURRENT_PED = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/family_relationships.txt"

# trim out CNVs, and variants which are present because they are compound hets
# with a CNV
exclude_cnv_based_variants <- function(variants) {
    
    cnvs = variants[grepl("DUP|DEL", variants$trio_genotype), c("proband", "chrom", "gene")]
    if (nrow(cnvs) == 0) { return(variants) }
    
    remove = cnvs[0, ]
    for (row_num in 1:nrow(cnvs)) {
        cnv = cnvs[row_num, ]
        cnv_genes = gsub(",", "|", cnv$gene)
        
        # get the variants for the same proband and on the same chrom
        matching = variants[variants$proband == cnv$proband & variants$chrom == cnv$chrom, ]
        
        # convert variant genes to vectors, and find which variants match any CNV gene
        genes = strsplit(matching$gene, ",")
        idx = sapply(genes, function(x) any(grepl(cnv_genes, x)))
        
        remove = rbind(remove, matching[idx, ])
    }
    
    # and drop the CNV related variants from the dataset
    variants = variants[!row.names(variants) %in% row.names(remove), ]
    
    return(variants)
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

#' run clinical filtering on a set of variants, debugging each time, so that can
#' identify the reason why each individual variant failed the clinical filtering
#'
#' @param variant data frame or list for single variant, providing proband,
#'     chrom and position info
#' @param ped_path path to a pedigree file that contains family information for
#'     the proband and the probands parents. Note that the ped file can contain
#'     info for other families, since we construct a ped file specifically for
#'     the proband of interest within this function, using the info from the
#'     ped_path argument.
#'
#' @examples
#' check_variant_fail(list(proband="DDDP100099", chrom=6, position=152712420),
#'    "~/samples.ped")
check_variant_fail <- function(variant, ped_path) {
    
    proband_id = variant[["proband"]]
    new_path = tempfile(fileext=".ped")
    construct_pedigree_file(proband_id, ped_path, new_path)
    
    pos = variant[["position"]]
    chrom = variant[["chrom"]]
    cat(c(proband_id, chrom, pos, "\t"))
    
    command = "python3"
    args = c("/nfs/users/nfs_j/jm33/apps/clinical-filter/src/main/python/clinical_filter.py",
    "--ped",  new_path,
    "--output", "test.txt",
    "--syndrome-regions", "/lustre/scratch113/projects/ddd/resources/decipher_syndrome_list_20140428.txt",
    "--known-genes", "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2014-11-04/DDG2P/dd_genes_for_clinical_filter",
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
    cat(c(reason, "\n"))
    file.remove(new_path)
    
    return(reason)
}

# aggregate the fail categories, so we can easily make tables
combine_fail_categories <- function(variants) {
    
    category = rep(NA, nrow(variants)
    
    # I've manually set all the failure categories. Most of the variants fail for
    # the MAF, polyphen prediction, and VQSLOD categories.
    category[grepl("LOW_VQSLOD", variants$reason)] = "failed LOW_VQSLOD"
    category[grepl("not smaller_than 0\\.01", variants$reason)] = "failed MAF"
    category[grepl("low MAF in non-biallelic", variants$reason)] = "failed MAF"
    category[grepl("variant not in VCF", variants$reason)] = "variant not in VCF"
    category[grepl("polyphen prediction", variants$reason)] = "polyphen changed"
    category[grepl("failed DENOVO-SNP\\/INDEL check", variants$reason)] = "not called by denovogear"
    category[grepl("gatk_LowQual", variants$reason)] = "failed gatk filter"
    category[grepl("gatk_StandardFilters", variants$reason)] = "failed gatk filter"
    category[grepl("TEAM29_FILTER", variants$reason)] = "failed denovogear filters"
    category[grepl("PP_DNM", variants$reason)] = "failed PP_DNM"
    category[grepl("failed CQ", variants$reason)] = "consequence changed"
    category[grepl("failed HGNC", variants$reason)] = "gene symbol not picked up"
    category[is.na(variants$reason)] = "likely missing compound pair"
    
    return(category)
}

main <- function() {
    diagnoses = read.xls(DIAGNOSES_PATH, sheet="Exome variants reviewed", stringsAsFactors=FALSE)
    diagnoses = subset(diagnoses, select=c("proband", "chrom", "position", "DECISION"))
    
    initial_path = file.path(HOME, "clinical_reporting.2014-12-10.txt")
    current_path = file.path(HOME, "clinical_reporting.2015-01-21.txt")
    
    initial = read.table(initial_path, sep="\t", stringsAsFactors=FALSE, header=TRUE)
    current = read.table(current_path, sep="\t", stringsAsFactors=FALSE, header=TRUE)
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
    
    # plot a venn diagram of the numbers of variants in each group
    v = venneuler(c(previous=nrow(initial_only), current=nrow(current_only), "previous&current"=(nrow(both)-(nrow(initial_only) + nrow(current_only)))))
    Cairo(file="variants_venn_diagram.pdf", type="pdf", height=10, width=10, units="cm")
    plot(v)
    dev.off()
    
    # check the reasons why each variant failed
    initial_only$reason = apply(initial_only, 1, check_variant_fail, CURRENT_PED)
    current_only$reason = apply(current_only, 1, check_variant_fail, ORIGINAL_PED)
    
    # categorize the fail reasons
    initial_only$category = combine_fail_categories(initial_only)
    current_only$category = combine_fail_categories(current_only)
    
    # get the numbers of variants in each fail category (in the inital_only
    # cross-tabulate by reporting decision)
    table(initial_only$category, initial_only$DECISION)
    table(current_only$category)
    
    # construct_pedigree_file("DDDP100149", "~/exome_reporting_2.ped", "temp.ped")
    # check_variant_fail(list(proband="DDDP100099", chrom=6, position=152712420), "~/samples.ped")
    
    # missing = read.table("~/problem_de_novos.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
    # missing$reason = apply(missing, 1, check_variant_fail, "~/exome_reporting_2.ped")
    # write.table(missing, file="~/problem_de_novos.txt", row.names=FALSE, sep="\t", quote=FALSE)
    
}


main()
