""" a helper script to process all the trios efficiently on the farm as a job
array. Splits trios in a PED file across multiple ped files in order to run
them in parallel.

Not recommended for use, as this is very scrappy code, and highly user-specific,
but it works if all the files are in the expected locations.
"""

import subprocess
import random
import os
import time
import glob
import argparse

USER_DIR = "~/"
APP_FOLDER = os.path.join(USER_DIR, "apps", "clinical-filter")

FILTER_CODE = os.path.join(APP_FOLDER, "src", "main", "python", "clinical_filter.py")

DATAFREEZE = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13"
KNOWN_GENES_PATH = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter"
ALT_IDS_PATH = os.path.join(DATAFREEZE, "person_sanger_decipher.txt")
PED_PATH = os.path.join(DATAFREEZE, "family_relationships.txt")
SYNDROMES_PATH = "/lustre/scratch113/projects/ddd/resources/decipher_syndrome_list_20140428.txt"
    
def get_options():
    """ gets the options from the command line
    """
    
    parser = argparse.ArgumentParser(description="Submit trio analysis to an \
        LSF cluster.")
    parser.add_argument("--ped", dest="ped_path", default=PED_PATH, \
        help="path of the ped file (default=construct from DDD datasets)")
    parser.add_argument("--log", dest="loglevel", \
        help="level of logging to use (default=all)")
    parser.add_argument("--ddg2p", dest="ddg2p_path", default=KNOWN_GENES_PATH, \
        help="optional path to the ddg2p file to use (default = current DDD DDG2P file)")
    parser.add_argument("--njobs", dest="n_jobs", default=100, \
        help="number of jobs you want to divide the run across")
    parser.add_argument("--all-genes", dest="all_genes", default=False, \
        action="store_true", help="Option to assess variants in all genes. \
        If unused, restricts variants to DDG2P genes.")
    parser.add_argument("--without-parents", default=False,
        action="store_true", help="whether to remove the parents and analyse \
        probands only")
    parser.add_argument("--use-singletons-with-parents", default=False,
        action="store_true", help="Some probands have parents, but genotypes \
        from one or both parents are not yet available. The genetic data for \
        these parents will eventually be available. Other probands lack one or \
        both parents (and they will never be available). The singleton probands \
        with parents can be distingushed in the ped file as probands where both \
        parental IDs are non-zero, but lines defining either parent do not \
        exist. We'll currently exclude the singletons with parents, as \
        filtering of their genetic variants would provide an abundance of \
        variants that would later be clarified with parental genotypes.")
    
    args = parser.parse_args()
    
    if args.without_parents:
        args.use_singletons_with_parents = True
    
    return args

def check_for_singletons_missing_parents(family):
    """ figure out whether a family has parents, but lack parental data
    
    Args:
        family: list of lines from ped for family members
    
    Returns:
        true/false for whether the family has a proband with defined parents,
        but the parents do not have VCFs available.
    """
    
    # make sure we are working on a copy of the family lines, so we don't
    # interfere with later processing
    family = family[:]
    
    # make sure all of the samples are split into lists
    for pos in range(len(family)):
        family[pos] = family[pos].strip().split("\t")
        
    # identify the members in the family
    children = []
    parents = []
    for sample in family:
        dad_id = sample[2]
        mom_id = sample[3]
        
        if dad_id != "0" or mom_id != "0":
            children.append(sample)
        else:
            parents.append(sample)
    
    # get the IDs for the mother and father (these are the same for all children
    # in the family).
    dad_id = children[0][2]
    mom_id = children[0][3]
    
    # probands who are missing one of their parents don't need to be excluded,
    # since their missing parent will never have genetic data available.
    if dad_id == "0" or mom_id == "0":
        return False
    
    return len(parents) != 2

def split_pedigree_file(tempname, ped_path, number_of_jobs, exclude_parents, use_singletons_with_parents):
    """ split the ped file into multiple smaller ped files
    
    Args:
        tempname: string for the output path
        ped_path: path to pedigree file
        number_of_jobs: how many computational jobs to split the families over.
            Note that due to how the families are striuctured (siblings etc), we
            might get more files than this
        exclude_parents: true/false for whether to exclude parents from the run.
        use_singletons_with_parents: true/false for whether to exclude probands
            who have parents define, but where one or both parents genotypes are
            not yet available.
    
    Returns:
        The number of files that the cohort has been split across (which will
        now be the number of jobs to run).
    """
    
    # create a dictionary of lines by family ID
    families = {}
    with open(ped_path, "r") as ped:
        for line in ped:
            split_line = line.split("\t")
            family_ID = split_line[0]
            
            # ignore header lines
            if line.startswith("family_id"):
                continue
            
            # if we want to exclude the parents, don't include lines with "0" for
            # the parental ID, these are parental lines.
            if exclude_parents and split_line[2] == "0":
                continue
            
            if family_ID not in families:
                families[family_ID] = []
            
            families[family_ID].append(line)
    
    if not use_singletons_with_parents:
        for family_id in sorted(families):
            family = families[family_id]
            if check_for_singletons_missing_parents(family):
                del families[family_id]
    
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
            output_file = open("{0}.{1}.txt".format(tempname, file_counter), "w")
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
    
    Returns:
        random 8 character string
    """
    
    # set up a random string to associate with the run
    hash_string = "{0:8x}".format(random.getrandbits(32))
    hash_string = hash_string.strip()
    
    # don't allow the random strings to be equivalent to a number, since
    # the LSF cluster interprets those differently from letter-containing
    # strings
    while is_number(hash_string):
        hash_string = "{0:8x}".format(random.getrandbits(32))
        hash_string = hash_string.strip()
    
    return hash_string

def submit_bsub_job(command, job_id=None, dependent_id=None, memory=None, requeue_code=None, logfile=None):
    """ construct a bsub job submission command
    
    Args:
        command: list of strings that form a unix command
        job_id: string for job ID for submission
        dependent_id: job ID, or list of job IDs which the current command needs
            to have finished before the current command will start. Note that
            the list can be empty, in which case there are no dependencies.
        memory: minimum memory requirements (in megabytes)
        requeue_code: exit codes under which it is acceptable to requeue the
            job, or None
        logfile: path for bjob output
    
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

def run_array(hash_string, n_jobs, temp_name, genes_path, all_genes, log_options):
    """ runs clinical filtering, split across multiple compute jobs
    
    Args:
        hash_string: random hash string as identifier for the jobs
        n_jobs: number of compute jobs to run (same as the number of ped files
            to run on, which we have created earlier)
        temp_name: prefix of the ped file basenames
        genes_path: path to the known genes file, or None
        all_genes: True/False whether to test all geens, or restrict variants to
            the known genes.
        log_options: logging level to use
    """
    
    # set up run parameters
    job_id = "{0}[1-{1}]".format(hash_string, n_jobs)
    bjob_output_name = "{0}.bjob_output".format(temp_name)
    output_name = "{0}.output".format(temp_name)
    
    command = ["python3", FILTER_CODE, \
        "--ped", "{0}.\$LSB_JOBINDEX\.txt".format(temp_name), \
        "--output", "{0}.\$LSB_JOBINDEX\.txt".format(output_name), \
        "--syndrome-regions", SYNDROMES_PATH, \
        "--alternate-ids", ALT_IDS_PATH, \
        "--pp-dnm-threshold", "0"] + log_options
    
    # sometimes we don't want to restrict to the DDG2P genes, then all_genes
    # would be False and variants would be assessed in every gene.
    if not all_genes:
        command += ["--known-genes", genes_path]
    
    submit_bsub_job(command, job_id=job_id, logfile=bjob_output_name + ".%I.txt")

def run_cleanup(hash_string):
    """ runs a lsf job to clean up the files
    """
    
    output_name = "tmp_ped.{0}.output".format(hash_string)
    
    # merge the array output after the array finishes
    merge_id = "merge1_{0}".format(hash_string)
    command = ["head", "-n", "1", output_name + ".1.txt", \
        ">", "clinical_reporting.txt", ";", \
        "tail", "-q", "-n", "+2", output_name + "*", \
        ">>", "clinical_reporting.txt"]
    submit_bsub_job(command, job_id=merge_id, dependent_id=hash_string, \
        logfile="tmp_ped.{0}*bjob_output.var_merge.txt".format(hash_string))
    
    # merge the log files after the array finishes
    log_merge_id = "merge2_{0}".format(hash_string)
    command = ["cat", "tmp_ped.{0}*.log".format(hash_string), ">", "clinical_reporting.log"]
    submit_bsub_job(command, job_id=log_merge_id, dependent_id=hash_string, \
        logfile="tmp_ped.{0}*bjob_output.log_merge.txt".format(hash_string))
    
    # remove the temporary files
    command = ["rm", "tmp_ped.{0}*".format(hash_string)]
    submit_bsub_job(command, job_id="cleanup", dependent_id=[merge_id, log_merge_id], \
        logfile="clinical_reporting.cleanup.txt")

def main():
    
    opts = get_options()
    
    log_options = []
    if opts.loglevel != None:
        log_options = ["--log", opts.loglevel]
    
    # set up a random string to associate with the run
    hash_string = get_random_string()
    
    temp_name = "tmp_ped.{0}".format(hash_string)
    
    tidy_directory_before_start()
    trio_counter = split_pedigree_file(temp_name, opts.ped_path, opts.n_jobs, opts.without_parents, opts.use_singletons_with_parents)
    run_array(hash_string, trio_counter, temp_name, opts.ddg2p_path, opts.all_genes, log_options)
    run_cleanup(hash_string)


if __name__ == "__main__":
    main()
