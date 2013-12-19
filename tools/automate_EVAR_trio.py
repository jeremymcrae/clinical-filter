'''
Given direcotry of Trio families (each subdir has 3 VCF files), this script
will automate EVAR_trio_0.3 and output the text results in the same folder.

sa9@sanger.ac.uk
2 July 2012
'''

import os, sys, fnmatch
import re

#target = "/Users/Macia/Desktop/DropBox/Dropbox/Dropbox/Team29/MyProjects/CHD/AVSD/Leuven/8_trios_23MAy12/vcf/merged_EPS"
target = sys.argv[1] # "/nfs/users/nfs_s/sa9/chd/Newcastle/VCF/3.uk10kTwin_ESP6500_annotated/"
#target = "/Users/Macia/Desktop/Sanger/TOF_vcfs/2011-11-03/raw_merged/Ready"

EVAR_trio = "../EVAR_trio_0.3.py"
filter= "../filters.txt"
hierarchy= "../hierarchy.txt"
output= "../output.txt"
tags= "../tags.txt"
children = []
reportType = sys.argv[2]

def loadChildrenList():
    global children
    TOF_children = '''TOF5135941
                TOF5135944
                TOF5135947
                TOF5135950
                TOF5135953
                TOF5135956
                TOF5135959
                TOF5135965
                TOF5135968
                TOF5135971
                TOF5135977
                TOF5135980
                TOF5135983
                TOF5135986
                TOF5135989
                TOF5135992
                TOF5135995
                TOF5135998
                TOF5136004
                TOF5136007
                TOF5136010
                TOF5136013
                TOF5136016
                TOF5136019
                TOF5136022
                TOF5136025
                TOF5136028
                TOF5165230
                TOF5165218
                TOF5165221'''
    children = [x.strip(" ") for x in TOF_children.split("\n")]
    pass

def isChild(f_name):
    for child in children:
        if re.search(child,f_name):
            return True
    return False

def main():
    loadChildrenList()
    for x in os.listdir(target):
        folder = os.path.join(target, x)
        #print folder
        lst = []
        if os.path.isdir(folder):
            for path, subdirs, files in os.walk(folder):
                for name in files:
                    if fnmatch.fnmatch(name, "*.vcf.gz"):
                        if len(lst) <= 3:
                            lst.append(os.path.join(path,name))
                        if len(lst) == 3:
                            parents = []
                            for sample in lst:
                                if isChild(sample):
                                    childSample = sample
                                else:
                                    parents.append(sample)
                            
                            if reportType == "FullReport":        
                                o_path = os.path.join(path, name.split(".")[0] + "_EVAR_FullReport.txt")
                                cmd = 'bsub -oo EVAR.log -M 2000000 -P cnpoly -R"select[mem>2000] rusage[mem=2000]" "python %s -c %s -f %s -m %s -l %s -r %s -o %s -t %s >  %s"' % (EVAR_trio,childSample,parents[0],parents[1], filter,hierarchy, output, tags, o_path)
                            elif reportType == "RReport":
                                o_path = os.path.join(path, name.split(".")[0] + "_EVAR_R.txt")
                                cmd = 'bsub -oo EVAR.log -M 2000000 -P cnpoly -R"select[mem>2000] rusage[mem=2000]" "python %s -c %s -f %s -m %s -l %s -r %s -o %s -t %s -s R >  %s"' % (EVAR_trio,childSample,parents[0],parents[1], filter,hierarchy, output, tags, o_path)
                            print cmd
                            print 50 * "-"
                            #os.system(cmd)
                            pass

main()
