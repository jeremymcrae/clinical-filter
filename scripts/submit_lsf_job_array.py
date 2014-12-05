""" a helper script to process all the trios efficiently on the farm as a job 
array. Splits trios in a PED file across multiple ped files in order to run 
them in parallel. 

Not recommended for use, as this is very scrappy code, and highly user-specific,
but it works if all the files are in the expected locations.
"""

import math
import subprocess
import random
import os
import shutil
import time
import glob
import argparse

home_folder = "/nfs/users/nfs_j/jm33/"
home_lustre = "/lustre/scratch113/teams/hurles/users/jm33/"
app_folder = os.path.join(home_folder, "apps", "clinical-filter")

filter_code = os.path.join(app_folder, "src", "main", "python", "clinical_filter.py")
filters = os.path.join(app_folder, "config", "filters.txt")
tag_names = os.path.join(app_folder, "config", "tags.txt")
deprecated_genes = os.path.join(app_folder, "config", "ddg2p_deprecated_hgnc.txt")

datafreeze = "/nfs/ddd0/Data/datafreeze/1133trios_20131218/"
known_genes = os.path.join(datafreeze, "DDG2P_with_genomic_coordinates_20131107_updated_TTN.tsv")
alternate_ids = os.path.join(datafreeze, "person_sanger_decipher.private.txt")
individuals_filename = os.path.join(datafreeze, "family_relationships.shared.txt")
working_vcfs_filename = os.path.join(datafreeze, "all_working_paths.private.txt")
syndrome_regions_filename = "/lustre/scratch113/projects/ddd/resources/decipher_syndrome_list_20140428.txt"
    
def get_options():
    """ gets the options from the command line
    """
    
    parser = argparse.ArgumentParser(description="Submit trio analysis to an \
        LSF cluster.")
    parser.add_argument('--ped', dest='ped_path', default=None, help='path of the ped file (default=construct from DDD datasets)')
    parser.add_argument('--log', dest='loglevel', help='level of logging to use (default=all)')
    parser.add_argument('--ddg2p', dest='ddg2p_path', default=known_genes, help='optional path to the ddg2p file to use (default = current DDD DDG2P file)')
    parser.add_argument('--njobs', dest='n_jobs', default=100, help='number of jobs you want to divide the run across')
    parser.add_argument('--all-genes', dest='all_genes', default=False, action="store_true", help='Option to assess variants in all genes. If unused, restricts variants to DDG2P genes.')
    
    args = parser.parse_args()
    
    return args

def make_ped(ped_filename):
    """ create a PED file for clinical filtering, using DDD datafreeze files
    
    The DDD datafreeze folder contains two files from which we can create a ped
    file for filtering, the first file contains the first five columns in 
    linkage pedigree format, the second file contains paths to VCF files for 
    all the individuals in the first file. Unfortunately, that second file 
    only contains the paths, which we have to search for the individual IDs.
    
    Args:
        ped_filename: path to write ped file output to.
    """
    
    individuals = open(individuals_filename, "r")
    working_vcfs = open(working_vcfs_filename, "r")
    
    # get all the paths to the VCF files
    vcfs = {}
    for line in working_vcfs:
        if line.strip() == "":
            continue
        split_line = line.strip().split("vcfs/")[1]
        individual_ID = split_line.split(".")[0]
        vcfs[individual_ID] = line
    
    output = open(ped_filename, "w")
    
    # for each individual in the proto ped file, get the VCF path, then write a new ped line
    for line in individuals:
        split_line = line.strip().split()
        individual_ID = split_line[1]
        if individual_ID == "individual_id":
            continue
        output.write(line.strip() + "\t" + vcfs[individual_ID])
    
    output.close()

def split_pedigree_file(tempname, ped_file, number_of_jobs):
    """ split the ped file into multiple smaller ped files
    """
    
    ped = open(ped_file, "r")
    
    # create a dictionary of lines by family ID
    families = {}
    for line in ped:
        family_ID = line[0:20].split("\t")[0]
        if family_ID not in families:
            families[family_ID] = []
        families[family_ID].append(line)
    
    # figure out how many families to include per file, in order to make the correct number of jobs
    max_trios_in_ped_file = float(len(families))/float(number_of_jobs)
    
    file_counter = 1
    family_counter = 0
    for family_ID in sorted(families):
        lines = families[family_ID]
        family_counter += 1
        if (file_counter * family_counter) > (file_counter * max_trios_in_ped_file):
            file_counter += 1
            family_counter = 1
        if family_counter == 1:
            output_file = open(tempname + str(file_counter) + ".txt", "w")
        output_file.writelines(lines)
        
    return file_counter

def tidy_directory_before_start():
    """ clean the folder of old job array files before we start our current run
    """
    
    # clean up the directory before we start
    files = glob.glob("tmp_ped*")
    for filename in files:
        os.remove(filename)
    
def write_sh_file(hash_string, command):
    """
    """
    
    # the command cannot be run directly from subprocess, as the quote marks are not processed 
    # correctly. Write it to a file, then check that directly.
    random_filename = hash_string + ".sh"
    random_file = open(random_filename, 'w')
    random_file.write(" ".join(command))
    random_file.close()
    
    return random_filename

def remove_sh_file(filename):
    """ cleans up the sh file
    """
    
    # wait for the script to be run and the job to be picked up, then remove the script file
    time.sleep(1)
    os.remove(filename)
    
    if os.path.exists("clinical_reporting.log"):
        os.remove("clinical_reporting.log")

def run_array(hash_string, trio_counter, temp_name, output_name, known_genes_path, all_genes, log_options):
    """ sets up a lsf job array
    """
    
    # set up run parameters
    job_name = hash_string + '[1-' + str(trio_counter)
    job_array_params = '-J "' + job_name + ']"'
    
    bjob_output_name = temp_name + "bjob_output"
    
    # copy the alternate IDs to a faster filesystem
    fast_alternate_ids = os.path.join(home_lustre, os.path.basename(alternate_ids))
    shutil.copyfile(alternate_ids, fast_alternate_ids)
    
    # copy the syndrome regions to potentiall faster system
    fast_syndrome_regions = os.path.join(home_lustre, os.path.basename(syndrome_regions_filename))
    shutil.copyfile(syndrome_regions_filename, fast_syndrome_regions)
    
    # command = ["bsub", job_array_params, "-o", bjob_output_name + ".%I.txt", "python3", filter_code, "--ped", temp_name + "\$LSB_JOBINDEX\.txt", "--filter", filters, "--tags", tag_names, "--alternate-ids", fast_alternate_ids, "--output", output_name + "\$LSB_JOBINDEX\.txt", "--syndrome-regions", syndrome_regions_filename] + log_options
    command = ["bsub", job_array_params, "-o", bjob_output_name + ".%I.txt", "python3", filter_code, "--ped", temp_name + "\$LSB_JOBINDEX\.txt", "--filter", filters, "--tags", tag_names, "--output", output_name + "\$LSB_JOBINDEX\.txt", "--syndrome-regions", fast_syndrome_regions, "--deprecated-genes", deprecated_genes] + log_options
    
    # sometimes we don't want to restrict to the DDG2P genes, then all_genes
    # would be False and variants would be assessed in every gene.
    if not all_genes:
        # copy the known genes to a faster filesystem
        genes_path = os.path.join(home_lustre, os.path.basename(known_genes_path))
        shutil.copyfile(known_genes_path, genes_path)
        command += ["--known-genes", genes_path]
    
    sh_file = write_sh_file(hash_string, command)
    subprocess.Popen(["sh", sh_file])
    remove_sh_file(sh_file)

def run_cleanup(hash_string, output_name, temp_name):
    """ runs a lsf job to clean up the files
    """
    
    # merge the array output after the array finishes
    command = ['bsub -J "merge1_' + hash_string + '" -w "done(' + hash_string + ')"', "-o", temp_name + "bjob_output.var_merge.txt", "bash", "-c", "\"head", "-n", "1", output_name + "1.txt", ">", "clinical_reporting.txt", "; tail", "-q", "-n", "+2", output_name + "*", ">>", "clinical_reporting.txt\""]
    sh_file = write_sh_file(hash_string, command)
    subprocess.Popen(["sh", sh_file])
    remove_sh_file(sh_file)
    
    # merge the log files after the array finishes
    command = ['bsub -J "merge2_' + hash_string + '" -w "done(' + hash_string + ')"', "-o", temp_name + "bjob_output.log_merge.txt", "bash", "-c", "\"cat", temp_name + "*.log", ">", "clinical_reporting.log\""]
    sh_file = write_sh_file(hash_string, command)
    subprocess.Popen(["sh", sh_file])
    remove_sh_file(sh_file)
    
    # submit a cleanup job to the cluster
    command = ['bsub -J "cleanup" -w "merge2_' + hash_string + '"', "-o", "clinical_reporting.cleanup.txt", "bash", "-c", "\"rm", temp_name + "*\""]
    sh_file = write_sh_file(hash_string, command)
    subprocess.Popen(["sh", sh_file])
    remove_sh_file(sh_file)

def is_number(string):
    """ check whether a string can be converted to a number
    """
    
    try:
        number = float(string)
    except ValueError:
        return False
    
    return True

def main():
    
    opts = get_options()
    
    if opts.loglevel != None:
        loglevel = opts.loglevel
        log_options = ["--log", loglevel]
    else:
        log_options = []
    
    if opts.ped_path is None:
        ped_path = os.path.join(home_folder, "exome_reporting.ped")
        make_ped(ped_path)
    else:
        ped_path = opts.ped_path
    
    # set up a random string to associate with the run
    hash_string = "%8x" % random.getrandbits(32)
    hash_string = hash_string.strip()
    
    while is_number(hash_string):
        hash_string = "%8x" % random.getrandbits(32)
        hash_string = hash_string.strip()
    
    temp_name = "tmp_ped." + hash_string + "."
    output_name = temp_name + "output."
    
    tidy_directory_before_start()
    trio_counter = split_pedigree_file(temp_name, ped_path, opts.n_jobs)
    run_array(hash_string, trio_counter, temp_name, output_name, opts.ddg2p_path, opts.all_genes, log_options)
    run_cleanup(hash_string, output_name, temp_name)


if __name__ == "__main__":
    main()



