"""" a helper script to process all the trios efficiently on the farm as a job
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

from clinicalfilter.ped import load_families

USER_DIR = os.path.expanduser('~')
APP_FOLDER = os.path.join(USER_DIR, "apps", "clinical-filter")

FILTER_CODE = os.path.join(APP_FOLDER, "bin", "clinical_filter.py")

DATAFREEZE = "/nfs/ddd0/ddd_scratch/final_8k_resources_candidate"
KNOWN_GENES = "/nfs/ddd0/ddd_scratch/ddg2p_prep/candidate_clinical_filter_ddg2p_file.txt"
PED_PATH = os.path.join(DATAFREEZE, "family_relationships.ped")
SYNDROMES_PATH = "/lustre/scratch113/projects/ddd/resources/decipher_syndrome_list_20140428.txt"
LAST_BASE_PATH = "/lustre/scratch113/projects/ddd/users/jm33/last_base_sites_G.json"
    
def get_options():
    """ gets the options from the command line
    """
    
    parser = argparse.ArgumentParser(description="Submit trio analysis to an \
        LSF cluster.")
    parser.add_argument("--ped", default=PED_PATH,
        help="path of the ped file (default=construct from DDD datasets)")
    parser.add_argument("--log", dest="loglevel",
        help="level of logging to use (default=all)")
    parser.add_argument("--ddg2p", default=KNOWN_GENES_PATH,
        help="optional path to the ddg2p file to use (default = current DDD DDG2P file)")
    parser.add_argument("--njobs", default=100,
        help="number of jobs you want to divide the run across")
    parser.add_argument("--all-genes", dest="all_genes", default=False,
        action="store_true", help="Option to assess variants in all genes. \
        If unused, restricts variants to DDG2P genes.")
    parser.add_argument("--without-parents", default=False, action="store_true",
        help="whether to remove the parents and analyse probands only")
    parser.add_argument("--use-singletons-with-parents", default=False,
        action="store_true", help="Some probands have parents, but genotypes"
        "from one or both parents are not yet available. The genetic data for"
        "these parents will eventually be available. Other probands lack one or"
        "both parents (and they will never be available). The singleton probands"
        "with parents can be distingushed in the ped file as probands where both"
        "parental IDs are non-zero, but lines defining either parent do not"
        "exist. We'll currently exclude the singletons with parents, as"
        "filtering of their genetic variants would provide an abundance of"
        "variants that would later be clarified with parental genotypes.")
    parser.add_argument("--ignore-lof-tweak", default=False, action="store_true",
        help="whether to use the last base of exon rule.")
    
    args = parser.parse_args()
    
    if args.without_parents:
        args.use_singletons_with_parents = True
    
    return args

def is_singleton_without_parents(family):
    """ find probands who have parents recruited and currently lack parental data
    
    Args:
        family: Family object
    
    Returns:
        true/false for whether the family has a proband with defined parents,
        but the parents do not have VCFs available.
    """
    
    child = family.children[0]
    return child.mom_id != '0' and child.dad_id != '0' and \
        (family.mother is None or family.father is None)

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
    
    families = load_families(ped_path)
    
    if not use_singletons_with_parents:
        families = [ x for x in families if not is_singleton_without_parents(x) ]
    
    # figure out how many families to include per file, in order to make the
    # correct number of jobs
    max_families = float(len(families))/float(number_of_jobs)
    
    files_n = 1
    families_n = 0
    for family in families:
        families_n += 1
        
        if families_n > max_families:
            files_n += 1
            families_n = 1
        
        if families_n == 1:
            output_file = open("{0}.{1}.txt".format(tempname, files_n), "w")
        
        for person in family:
            if person is None:
                continue
            
            line = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(person.family_id,
                person.get_id(), person.dad_id, person.mom_id,
                person.get_gender(), person.get_affected_status(),
                person.get_path())
            output_file.write(line)
    
    return files_n

def clean_folder(string):
    ''' clean the folder of old job array files'''
    
    for path in glob.glob(string + "*"):
        os.remove(path)

def get_random_string():
    """ make a random string, which we can use for bsub job IDs, so that
    different jobs do not have the same job IDs.
    
    Returns:
        random 8 character string
    """
    
    def is_number(string):
        try:
            number = float(string)
            return True
        except ValueError:
            return False
    
    # don't allow the random strings to be equivalent to a number, since
    # the LSF cluster interprets those differently from letter-containing
    # strings
    job_id = None
    while job_id is None or is_number(job_id):
        job_id = "{0:x}".format(random.getrandbits(32))
        job_id = job_id.strip()
    
    return job_id

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

def run_array(hash_string, n_jobs, temp_name, genes_path, all_genes, ignore_lof_tweak, log_options):
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
    
    command = ["python3", FILTER_CODE,
        "--ped", "{0}.\$LSB_JOBINDEX\.txt".format(temp_name),
        "--output", "{0}.\$LSB_JOBINDEX\.txt".format(output_name),
        "--syndrome-regions", SYNDROMES_PATH,
        "--pp-dnm-threshold", "0"] + log_options
    
    # sometimes we don't want to restrict to the DDG2P genes, then all_genes
    # would be False and variants would be assessed in every gene.
    if not all_genes:
        command += ["--known-genes", genes_path]
    
    if not ignore_lof_tweak:
        command += ["--lof-sites", LAST_BASE_PATH]
    
    submit_bsub_job(command, job_id=job_id, logfile=bjob_output_name + ".%I.txt", memory=150)

def run_cleanup(hash_string):
    """ runs a lsf job to clean up the files
    """
    
    date = time.strftime("%Y-%m-%d", time.localtime())
    output_name = "tmp_ped.{0}.output".format(hash_string)
    merged_output = "clinical_reporting.{}.txt".format(date)
    merged_log = "clinical_reporting.{}.log".format(date)
    
    # merge the array output after the array finishes
    merge_id = "merge1_{0}".format(hash_string)
    command = ["head", "-n", "1", output_name + ".1.txt",
        ">", merged_output, ";",
        "tail", "-q", "-n", "+2", output_name + "*",
        ">>", merged_output]
    submit_bsub_job(command, job_id=merge_id, dependent_id=hash_string,
        logfile="tmp_ped.{0}*bjob_output.var_merge.txt".format(hash_string))
    
    # merge the log files after the array finishes
    log_merge_id = "merge2_{0}".format(hash_string)
    command = ["cat", "tmp_ped.{0}*.log".format(hash_string), ">", merged_log]
    submit_bsub_job(command, job_id=log_merge_id, dependent_id=hash_string,
        logfile="tmp_ped.{0}*bjob_output.log_merge.txt".format(hash_string))
    
    # remove the temporary files
    command = ["rm", "tmp_ped.{0}*".format(hash_string)]
    submit_bsub_job(command, job_id="cleanup", dependent_id=[merge_id, log_merge_id],
        logfile="clinical_reporting.cleanup.txt")

def main():
    
    args = get_options()
    
    log_options = []
    if args.loglevel != None:
        log_options = ["--log", args.loglevel]
    
    # set up a random string to associate with the run
    hash_string = get_random_string()
    
    temp_name = "tmp_ped.{0}".format(hash_string)
    
    clean_folder('tmp_ped')
    trio_counter = split_pedigree_file(temp_name, args.ped, args.njobs,
        args.without_parents, args.use_singletons_with_parents)
    run_array(hash_string, trio_counter, temp_name, args.ddg2p, \
        args.all_genes, args.ignore_lof_tweak, log_options)
    run_cleanup(hash_string)


if __name__ == "__main__":
    main()
