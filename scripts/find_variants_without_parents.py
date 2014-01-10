""" script to make a pedigree file for individuals without their parents, and
then submit a cluster array job filtering on the PED file.
"""


import os
import sys
import optparse
import subprocess
import random
import glob

user_folder = "/nfs/users/nfs_j/jm33/"
app_folder = os.path.join(home_folder, "apps", "clinical-filter", "src", "main", "python")
submit_script = os.path.join(app_folder, "scripts", "submit_lsf_job_array.py")
ped_file = os.path.join(user_folder, "exome_reporting.ped")

def get_options():
    """ gets the options from the command line
    """
    
    parser = optparse.OptionParser()
    parser.add_option('--log', dest='loglevel', default="debug", help='level of logging to use, choose from: debug, info, warning, error or critical')
    
    (opts, args) = parser.parse_args()
    
    return opts

def load_ped(ped_path):
    """ loads a pedigree file and excludes any parents
    """
    
    ped = open(ped_path, "r")
    ped_read = ped.readlines()
    ped.close()
    
    # find the maternal and paternal IDs for the individual
    new_ped_lines = []
    for line in ped_read:
        split_line = line.split()
        individual_ID = split_line[1]
        paternal_ID = split_line[2]
        maternal_ID = split_line[3]
        
        if paternal_ID != "0" and maternal_ID != "0":
            new_ped_lines.append(line)
    
    return new_ped_lines

def main():
    """ analyse a single proband using bjobs
    """
    
    options = get_options()
    logging_option = ["--log", options.loglevel]
    
    new_ped = load_ped(ped_file)
    
    # remove the temp files from the previous run
    tmp_name = "tmp_run."
    files = glob.glob(tmp_name + "*")
    for filename in files:
        os.remove(filename)
    
    # write the temp ped file for the family to a file, but make sure it doesn't overwrite anything
    random_filename = tmp_name + str(random.randint(100000, 999999)) + ".ped"
    while os.path.exists(random_filename):
        random_filename = tmp_name + str(random.randint(100000, 999999)) + ".ped"
    
    random_file = open(random_filename, 'w')
    random_file.writelines(new_ped)
    random_file.close()
    
    # now set up the command for analysing the given pedigree file
    filter_command = ["python", submit_script, "--ped", random_filename] + logging_option
    
    subprocess.call(filter_command)

if __name__ == "__main__":
    main()





