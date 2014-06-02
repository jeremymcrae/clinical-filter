#/bin/sh

proband_id=$1
position=$2

# define the files to examine
datafreeze_dir=/nfs/ddd0/Data/datafreeze/1133trios_20131218
family_file="$datafreeze_dir"/family_relationships.shared.txt
paths_file="$datafreeze_dir"/all_working_paths.private.txt

# find the IDs of the parents
father_id=`grep "$proband_id" "$family_file" | cut -f3`
mother_id=`grep "$proband_id" "$family_file" | cut -f4`

function el_index {
    cnt=0; for el in "${ar[@]}"; do
        [[ $el == "$1" ]] && echo $cnt && break
        ((++cnt))
    done
}

for sample_id in "$proband_id" "$mother_id" "$father_id"; 
do
    # get the paths to the vcf
    sample_path=`grep "$sample_id" "$paths_file"`
    
    # find the format keys and values for the variant
    vcf_line=`zgrep "$position" "$sample_path"`
    # echo "$vcf_line"
    depths=`echo "$vcf_line" | cut -f9-10`
    ref=`echo "$vcf_line" | cut -f4`
    alt=`echo "$vcf_line" | cut -f5`
    
    # check if the var is an indel
    indel="SNV"
    if (( ${#ref} > 1 || ${#alt} > 1)) ; then
        indel="INDEL"
    fi
    
    # get an array of format keys, find the read depth and allele depth indices
    format=`echo $depths | cut -d ' '  -f1`
    ar=(${format//:/ })
    read_depth_index=`el_index "gatk_DP"`
    allele_depth_index=`el_index "AD"`
    
    # and get the read depth and allele depth values from the indices
    format=`echo $depths | cut -d ' '  -f2`
    ar=(${format//:/ })
    echo -e "$sample_id\t${position}\t${indel}\t${ar["$read_depth_index"]}\t${ar["$allele_depth_index"]}"
done

