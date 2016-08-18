""" runs clinical filtering analysis on a single proband, just starting from the
proband ID. A helper script used for testing purposes.

Not recommended for use, as this is very scrappy code, and highly user-specific,
but it works if all the files are in the expected locations.
"""

import os
import sys
import argparse
import subprocess
import random
import glob


home_folder = "/nfs/users/nfs_j/jm33/"
app_folder = os.path.join(home_folder, "apps", "clinical-filter")

filter_code = os.path.join(app_folder, "bin", "clinical_filter.py")

datafreeze = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2015-04-13"
known_genes = "/lustre/scratch113/projects/ddd/resources/ddd_data_releases/2015-04-13/DDG2P/dd_genes_for_clinical_filter"
ped_file = os.path.join(datafreeze, "family_relationships.txt")
alternate_ids = os.path.join(datafreeze, "person_sanger_decipher.txt")
syndrome_regions_filename = "/lustre/scratch113/projects/ddd/resources/decipher_syndrome_list_20140428.txt"
LAST_BASE_PATH = "/lustre/scratch113/projects/ddd/users/jm33/last_base_sites_G.json"

def get_options():
    """ gets the options from the command line
    """
    
    parser = argparse.ArgumentParser(description="Submit analysis job for single individual")
    parser.add_argument('-i', '--individual', dest='proband_ID', required=True, help='ID of proband to be analysed')
    parser.add_argument('--log', dest='loglevel', default="debug", help='level of logging to use, choose from: debug, info, warning, error or critical')
    parser.add_argument('--all-genes', dest='all_genes', default=False, action="store_true", help='Option to assess variants in all genes. If unused, restricts variants to DDG2P genes.')
    parser.add_argument('--debug-chrom', dest='debug_chrom', help='chromosome of variant to debug.')
    parser.add_argument('--debug-pos', dest='debug_pos', help='position of variant to debug.')
    parser.add_argument('--without-parents', default=False, action="store_true",
        help='whether to remove the parents for a proband only-analysis.')
    
    args = parser.parse_args()
    
    return args

def load_ped(ped_path, proband_ID, exclude_parents):
    """ loads the pedigree details for a prband
    
    Args:
        ped_path: path to pedigree file for cohort
        proband_ID: individual Id for proband of interest
        exclude_parents: whether to exclude the parents of the proband
    """
    
    with open(ped_path, "r") as ped:
        ped_read = ped.readlines()
    
    # find the maternal and paternal IDs for the individual
    maternal_ID = None
    paternal_ID = None
    for line in ped_read:
        line = line.split()
        individual_ID = line[1]
        
        if individual_ID == proband_ID:
            family_ID = line[0]
            paternal_ID = line[2]
            maternal_ID = line[3]
            break
    
    if maternal_ID is None:
        sys.exit("Could not find sample ID in family relationships file.")
    
    # extract all the lines for the individual
    new_ped_lines = []
    for line in ped_read:
        split_line = line.strip().split()
        individual_ID = split_line[1]
        
        if individual_ID == proband_ID:
            new_ped_lines.append(line)
        elif not exclude_parents and (individual_ID == paternal_ID or individual_ID == maternal_ID):
            new_ped_lines.append(line)
    
    return new_ped_lines


def main():
    """ analyse a single proband using bjobs
    """
    
    options = get_options()
    proband_ID = options.proband_ID
    logging_option = ["--log", options.loglevel]
    
    new_ped = load_ped(ped_file, proband_ID, options.without_parents)
    
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
    filter_command = ["python3", filter_code, \
        "--ped", random_filename, \
        "--alternate-ids", alternate_ids, \
        "--output", random_filename + ".output.txt", \
        "--export-vcf", os.getcwd(), \
        "--syndrome-regions", syndrome_regions_filename,
        "--lof-sites", LAST_BASE_PATH] + logging_option
    
    if not options.all_genes:
        filter_command += ["--known-genes", known_genes]
    
    if options.debug_chrom is not None:
        filter_command += ["--debug-chrom", options.debug_chrom, "--debug-pos", options.debug_pos]
        
    full_command = " ".join(bjobs_preamble + filter_command)
    
    subprocess.call(bjobs_preamble + filter_command)

if __name__ == "__main__":
    main()
