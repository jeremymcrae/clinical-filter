""" runs clinical filtering analysis on a single proband, just starting from the
proband ID. A helper script used for testing purposes.

Not recommended for use, as this is very scrappy code, and highly user-specific,
but it works if all the files are in the expected locations.
"""

import os
import argparse
import subprocess

from clinicalfilter.ped import load_families

HOME = os.path.expanduser('~')
SCRIPT = os.path.join(HOME, 'apps', "clinical-filter", "bin", "clinical_filter.py")

KNOWN_GENES = "/lustre/scratch115/projects/ddd/users/jm33/ddg2p.for_8K.with_hgnc_ids.txt"
FULL_PED = "/nfs/ddd0/Data/datafreeze/ddd_data_releases/2016-10-03/family_relationships.ped"
SYNDROMES = "/lustre/scratch113/projects/ddd/resources/decipher_syndrome_list_20140428.txt"
LAST_BASE_PATH = "/lustre/scratch113/projects/ddd/users/jm33/last_base_sites_G.json"

def get_options():
    """ gets the options from the command line
    """
    
    parser = argparse.ArgumentParser(description="Analyse a single individual")
    parser.add_argument('-i', '--individual', required=True,
        help='ID of proband to be analysed')
    parser.add_argument('--ped', default=FULL_PED, help='pedigree file to use')
    parser.add_argument("--known-genes", default=KNOWN_GENES,
        help="optional path to the known-genes file")
    parser.add_argument('--all-genes', default=False, action="store_true",
        help='Whether to check in all genes, regardless of known status.')
    parser.add_argument('--debug-chrom', help='chromosome of variant to debug.')
    parser.add_argument('--debug-pos', help='position of variant to debug.')
    
    return parser.parse_args()

def load_ped(ped_path, proband_id):
    """ loads the pedigree details for a prband
    
    Args:
        ped_path: path to pedigree file for cohort
        proband_ids: list of person_ids for probands of interest
    """
    
    families = load_families(ped_path)
    families = [ f for f in families for x in f.children if x.get_id() == proband_id ]
    family = families[0]
    
    to_line = lambda x: '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(x.family_id, x.get_id(),
        x.dad_id, x.mom_id, x.get_gender(), x.get_affected_status(), x.get_path())
    
    return [ to_line(x) for x in family if x is not None ]

def main():
    """ analyse a single proband using bjobs
    """
    
    args = get_options()
    
    ped = load_ped(args.ped, args.individual)
    
    # write the ped to a file
    path = "{}.ped".format(args.individual)
    with open(path, 'w') as handle:
        handle.writelines(ped)
    
    # now set up the command for analysing the given pedigree file
    command = ["python", SCRIPT, \
        "--ped", path, \
        "--output", "{}.output.txt".format(args.individual), \
        "--export-vcf", os.getcwd(), \
        "--syndrome-regions", SYNDROMES]
    
    if not args.all_genes:
        command += ["--known-genes", args.known_genes]
    
    if args.debug_chrom is not None:
        command += ["--debug-chrom", args.debug_chrom, "--debug-pos", args.debug_pos]
    
    subprocess.call(command)

if __name__ == "__main__":
    main()
