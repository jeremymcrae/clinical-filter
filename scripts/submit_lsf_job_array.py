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

datafreeze = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/"
known_genes = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2014-11-04/DDG2P/dd_genes_for_clinical_filter"
alternate_ids = os.path.join(datafreeze, "person_sanger_decipher.txt")
individuals_filename = os.path.join(datafreeze, "family_relationships.txt")
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
    individuals.readline()
    
    output = open(ped_filename, "w")
    for line in individuals:
        output.write(line)
    
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
            output_file = open("{0}{1}.txt".format(tempname, file_counter), "w")
        output_file.writelines(lines)
        
    return file_counter

def tidy_directory_before_start():
    """ clean the folder of old job array files before we start our current run
    """
    
    # clean up the directory before we start
    files = glob.glob("tmp_ped*")
    for filename in files:
        os.remove(filename)

def is_number(string):
    """ check whether a string can be converted to a number
    
    Args:
        string: value as a string, could be a number
        
    Returns:
        True/False for whether the value can be converted to a number
    """
    
    try:
        number = float(string)
    except ValueError:
        return False
    
    return True

def get_random_string():
    """ make a random string, which we can use for bsub job IDs, so that 
    different jobs do not have the same job IDs.
    """
    
    # set up a random string to associate with the run
    hash_string = "%8x" % random.getrandbits(32)
    hash_string = hash_string.strip()
    
    # done't allow the random strings to be equivalent to a number, since 
    # the LSF cluster interprets those differently from letter-containing 
    # strings
    while is_number(hash_string):
        hash_string = "%8x" % random.getrandbits(32)
        hash_string = hash_string.strip()
    
    return hash_string


def submit_bsub_job(command, job_id=None, dependent_id=None, memory=None, requeue_code=None, logfile=None):
    """ construct a bsub job submission command
    
    Args:
        command: list of strings that forma unix command
        job_id: string for job ID for submission
        dependent_id: job ID, or list of job IDs which the current command needs
            to have finished before the current command will start. Note that 
            the list can be empty, in which case there are no dependencies.
        memory: minimum memory requirements (in megabytes)
    
    Returns:
        nothing
    """
    
    if job_id is None:
        job_id = get_random_string()
    
    job = "-J \"{0}\"".format(job_id)
    
    mem = ""
    if memory is not None:
        mem = "-R 'select[mem>{0}] rusage[mem={0}]' -M {0}".format(memory)
    
    requeue = ""
    if requeue_code is not None:
        requeue = "-Q 'EXCLUDE({0})'".format(requeue_code)
    
    dependent = ""
    if dependent_id is not None:
        if type(dependent_id) == list:
            dependent_id = " && ".join(dependent_id)
        dependent = "-w '{0}'".format(dependent_id)
    
    log = "bjob_output.txt"
    if logfile is not None:
        log = logfile
    
    preamble = ["bsub", job, dependent, requeue, "-q", "normal", "-o", log, mem]
    command = ["bash", "-c", "\""] + command + ["\""]
    
    command = " ".join(preamble + command)
    subprocess.call(command, shell=True)

def run_array(hash_string, trio_counter, temp_name, output_name, known_genes_path, all_genes, log_options):
    """ sets up a lsf job array
    """
    
    # set up run parameters
    job_id = "{0}[1-{1}]".format(hash_string, trio_counter)
    bjob_output_name = temp_name + "bjob_output"
    
    # copy the alternate IDs to a faster filesystem
    fast_alternate_ids = os.path.join(home_lustre, os.path.basename(alternate_ids))
    shutil.copyfile(alternate_ids, fast_alternate_ids)
    
    # copy the syndrome regions to potentially faster filesystem
    fast_syndrome_regions = os.path.join(home_lustre, os.path.basename(syndrome_regions_filename))
    shutil.copyfile(syndrome_regions_filename, fast_syndrome_regions)
    
    command = ["python3", filter_code, \
        "--ped", "{0}\$LSB_JOBINDEX\.txt".format(temp_name), \
        "--output", "{0}\$LSB_JOBINDEX\.txt".format(output_name), \
        "--syndrome-regions", fast_syndrome_regions] + log_options
    
    # sometimes we don't want to restrict to the DDG2P genes, then all_genes
    # would be False and variants would be assessed in every gene.
    if not all_genes:
        # copy the known genes to a faster filesystem
        genes_path = os.path.join(home_lustre, os.path.basename(known_genes_path))
        shutil.copyfile(known_genes_path, genes_path)
        command += ["--known-genes", genes_path]
    
    submit_bsub_job(command, job_id=job_id, logfile=bjob_output_name + ".%I.txt")

def run_cleanup(hash_string, output_name, temp_name):
    """ runs a lsf job to clean up the files
    """
    
    # merge the array output after the array finishes
    merge_id = "merge1_{0}".format(hash_string)
    command = ["head", "-n", "1", output_name + "1.txt", ">", "clinical_reporting.txt", "; tail", "-q", "-n", "+2", output_name + "*", ">>", "clinical_reporting.txt"]
    submit_bsub_job(command, job_id=merge_id, dependent_id=hash_string, logfile=temp_name + "bjob_output.var_merge.txt")
    
    # merge the log files after the array finishes
    log_merge_id = "merge2_" + hash_string
    command = ["cat", temp_name + "*.log", ">", "clinical_reporting.log"]
    submit_bsub_job(command, job_id=log_merge_id, dependent_id=hash_string, logfile=temp_name + "bjob_output.log_merge.txt")
    
    # remove the temporary files
    command = ["rm", temp_name + "*"]
    submit_bsub_job(command, job_id="cleanup", dependent_id=[merge_id, log_merge_id], logfile="clinical_reporting.cleanup.txt")

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
    hash_string = get_random_string()
    
    temp_name = "tmp_ped.{0}.".format(hash_string)
    output_name = temp_name + "output."
    
    tidy_directory_before_start()
    trio_counter = split_pedigree_file(temp_name, ped_path, opts.n_jobs)
    run_array(hash_string, trio_counter, temp_name, output_name, opts.ddg2p_path, opts.all_genes, log_options)
    run_cleanup(hash_string, output_name, temp_name)


if __name__ == "__main__":
    main()



