""" script to make a pedigree file for individuals without their parents, and
then submit a cluster array job filtering on the PED file.
"""


import os
import sys
import argparse
import subprocess
import random
import glob

user_folder = "/nfs/users/nfs_j/jm33/"
app_folder = os.path.join(user_folder, "apps", "clinical-filter")
submit_script = os.path.join(app_folder, "scripts", "submit_lsf_job_array.py")
DATAFREEZE = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2014-11-04/"
PED_PATH = os.path.join(DATAFREEZE, "family_relationships.txt")

def get_options():
    """ gets the options from the command line
    """
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--log', dest='loglevel', default="debug", \
        help='level of logging to use, choose from: debug, info, warning, error or critical')
    parser.add_argument('--all-genes', dest='all_genes', default=False, \
        action="store_true", help='Option to assess variants in all genes. If unused, restricts variants to DDG2P genes.')
    
    args = parser.parse_args()
    
    return args

def load_ped(ped_path):
    """ loads a pedigree file and excludes any parents, leaving only probands
    """
    
    # write the temp ped file for the family to a file, but make sure it doesn't overwrite anything
    new_path = "tmp_run.{0}.ped".format(random.randint(100000, 999999))
    while os.path.exists(new_path):
        new_path = "tmp_run.{0}.ped".format(random.randint(100000, 999999))
    
    new_ped = open(new_path, 'w')
    
    # find the proband only lines, and write them to the new file
    with open(ped_path) as handle:
        for line in handle:
            split_line = line.split()
            individual_ID = split_line[1]
            paternal_ID = split_line[2]
            maternal_ID = split_line[3]
            
            if paternal_ID != "0" and maternal_ID != "0":
                new_ped.write(line)
    
    new_ped.close()
    
    return new_path

def tidy_directory_before_start():
    """ clean the folder of old job array files before we start our current run
    """
    
    # clean up the directory before we start
    files = glob.glob("tmp_run*")
    for filename in files:
        os.remove(filename)

def main():
    """ analyse a single proband using bjobs
    """
    
    options = get_options()
    logging_option = ["--log", options.loglevel]
    
    all_genes_option = []
    if options.all_genes:
        all_genes_option = ["--all-genes"]
    
    tidy_directory_before_start()
    new_ped = load_ped(PED_PATH)
    
    # now set up the command for analysing the given pedigree file
    filter_command = ["python", submit_script, "--ped", new_ped] + logging_option + all_genes_option
    
    subprocess.call(filter_command)

if __name__ == "__main__":
    main()
