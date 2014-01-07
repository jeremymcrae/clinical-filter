""" runs EVAR analysis on a single proband, just starting from the proband ID. Used for testing
purposes.
"""

import os
import sys
import optparse
import subprocess
import random
import glob


home_folder = "/nfs/users/nfs_j/jm33/"
app_folder = os.path.join(home_folder, "apps", "clinical-filter")

evar_code = os.path.join(app_folder, "clinical-filter.py")
ped_file = os.path.join(home_folder, "exome_reporting.ped")
filters = os.path.join(app_folder, "config", "filters.txt")
tag_names = os.path.join(app_folder, "config", "tags.txt")

datafreeze = "/nfs/ddd0/Data/datafreeze/1139trios_20131030/"
known_genes = os.path.join(datafreeze, "DDG2P_with_genomic_coordinates_20131107_updated_TTN.tsv")
alternate_ids = os.path.join(datafreeze, "person_sanger_decipher.private.txt")


def get_options():
    """ gets the options from the command line
    """
    
    parser = optparse.OptionParser()
    parser.add_option('-i', '--individual', dest='proband_ID', help='ID of proband to be analysed')
    parser.add_option('--log', dest='loglevel', default="debug", help='level of logging to use, choose from: debug, info, warning, error or critical')
    
    
    (opts, args) = parser.parse_args()
    
    return opts

def load_ped(ped_path, proband_ID):
    """ loads the pedigree details for a prband
    """
    
    ped = open(ped_path, "r")
    ped_read = ped.readlines()
    ped.close()
    
    # find the maternal and paternal IDs for the individual
    for line in ped_read:
        line = line.split()
        individual_ID = line[1]
        
        if individual_ID == proband_ID:
            family_ID = line[0]
            paternal_ID = line[2]
            maternal_ID = line[3]
            break
    
    # extract all the lines for the individual
    new_ped_lines = []
    for line in ped_read:
        split_line = line.strip().split()
        individual_ID = split_line[1]
        
        if individual_ID == proband_ID or individual_ID == paternal_ID or individual_ID == maternal_ID:
            new_ped_lines.append(line)
    
    return new_ped_lines


def main():
    """ analyse a single proband using bjobs
    """
    
    options = get_options()
    proband_ID = options.proband_ID
    logging_option = ["--log", options.loglevel]
    
    if proband_ID == None:
        sys.exit("you need to specify a proband ID using '-i' or '--individual'")
    
    new_ped = load_ped(ped_file, proband_ID)
    
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
    bjobs_preamble = ["bsub", "-q", "normal", "-o", random_filename + ".bjob_output.txt"]
    evar_command = ["python3", evar_code, "--ped", random_filename, "--filter", filters, "--tags", tag_names, "--known-genes", known_genes, "--alternate-ids", alternate_ids, "--output", random_filename + ".output.txt", "--export-vcf"] + logging_option
    full_command = " ".join(bjobs_preamble + evar_command)
    
    subprocess.call(bjobs_preamble + evar_command)

if __name__ == "__main__":
    main()







