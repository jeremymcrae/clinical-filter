'''
EVAR (Exome Variant Analysis Report) for single affected patiant

Given one VCF file for an affected patient, return  variants under AR , AR_dbHET
and AD.
step 1) create search space by using filters from filters.txt
step 2) create trees based on user defined hierarchy and boundaries in hierarchy.txt
step 3) print report
        a) basic information (data, input files, filters, hierarchies, etc )
        b) Counts of each group under AR, AR_dbHET and AD models
        c) Under each model (AR, AR_dbHET and AD) for each group , print its
           variants (only column defined by the user in output.txt)

Usage :

    python EVAR_single_0.1.py \
    -a single.vcf  \ 
    -f filters.txt \
    -o output.txt \
    -r hierarchy.txt \
    -t tags.txt \
    -s R \  # optional to print R-compatible format

Note:
The script accept vcf or vcf.gz files of a single or multiple samples as long
as the all cases or all controls not a mix.

version 0.1 (10 Jan 2012)

sa9@sanger.ac.uk
'''
import numpy as np
import os
import sys
import optparse

import user
import vcf
import parser


class Single(parser.Parser):
    def __init__(self,paths,weights_path, filters_path, output_path, tags_path, save):

        self.weights_path = weights_path
        self.filters_path = filters_path
        self.output_path = output_path
        self.paths = paths
        
        self.weights = user.parseWeights(weights_path)
        self.filters = user.parseFilters(filters_path)
        self.orders = user.parseOrders(output_path, weights_path)
        self.GN_tag, self.CQ_tag, self.MAF_tag, self.GT_tag = user.parseTags(tags_path)
        
        self.genes_dict = {}
        self.affectedTSVs = {}
        
        for path in paths:
            sample_name = os.path.basename(path)
            self.affectedTSVs[sample_name] = vcf.vcf2tsv(path)
            
        for sample_name in sorted(self.affectedTSVs.keys()):
            self.parse(self.affectedTSVs[sample_name], sample_name)
        
        self.createMatrix()
        self.get_var_under_models()
        
        if save == None:
            self.printReportInputInfo('\tSingle patient analysis [EVA Report]')
            self.printStatsTable()
            self.printResults()
        elif save == "R":
            self.save2R()
        elif save == "Excel":
            self.save2Excel()
    
    
    def createMatrix(self):
        print "#Creating genotypes and weights matrices"
        for gene in self.genes_dict.keys():
            sorted_pos = sorted([int(pos) for pos in self.genes_dict[gene]["POSs"].keys()])
            GTs = {}
            WTs = {}
            
            for member in sorted(self.affectedTSVs.keys()):
                if not GTs.has_key(member):
                    GTs[member] = [] # genotypes
                    WTs[member] = [] # weights

            for pos in sorted_pos:
                 for member in GTs.keys():
                    try:
                        record = self.genes_dict[gene]["POSs"][str(pos)][member]["raw"]
                        GT = vcf.translateGT(record[self.GT_tag])
                        WT = self.getWeights(record)
                        
                        GTs[member].append(GT)
                        WTs[member].append(WT)
                    except Exception as e:
                        GTs[member].append(0) # the member does not has variants in this pos i.e. hom ref
                        zeros = []
                        for key in self.weights.keys():
                            zeros.append(0.0)
                        WTs[member].append(zeros)

            self.genes_dict[gene]['GTs'] = np.array([GT_lst[1] for GT_lst in sorted(GTs.items())])
            self.genes_dict[gene]['WTs'] = np.array([WT_lst[1] for WT_lst in sorted(WTs.items())])
    
    def get_var_under_models(self):
        self.AR_results = self.construct_holder()
        self.AR_dbHET_results = self.construct_holder()
        self.AD_results = self.construct_holder()
        
        wt_length = len(self.orders["weights"])
        
        for member in sorted(self.affectedTSVs.keys()):
            header = self.affectedTSVs[member]['header']
            header.extend(['appears in','sample_name'])
            break
        
        for gene in self.genes_dict.keys():
            
            GTs = self.genes_dict[gene]['GTs']
            #print GTs
            # genotypes under compound /double HET needs speical treatment.
            # Copy them for further manipulation 
            GTs_db_HET = GTs
            
            WTs = self.genes_dict[gene]['WTs']
            sorted_pos = sorted([int(pos) for pos in self.genes_dict[gene]["POSs"].keys()])
            cln_idx = GTs.shape[1]
            
            for i in range(0,cln_idx):
                variant_GT = GTs[0,i]
                key = tuple(x for x in WTs[0,i]) # only include child weights , row 0                
                
                if variant_GT == 2:
                    record = self.extractUserFields(self.genes_dict[gene]['POSs'][str(sorted_pos[i])][member]['raw'], header)
                    
                    if not self.AR_results.has_key(key):
                        self.AR_results[key] = {gene:{str(sorted_pos[i]):record}}
                    else:
                        if not self.AR_results[key].has_key(gene):
                            self.AR_results[key][gene] = {str(sorted_pos[i]):record}
                        else:
                            self.AR_results[key][gene][str(sorted_pos[i])] = record
                    #print key, "\t".join(record)

                
                # for loci with at least two sibs with HET
                if variant_GT == 1:
                    record = self.extractUserFields(self.genes_dict[gene]['POSs'][str(sorted_pos[i])][member]['raw'], header)
                
                    if not self.AD_results.has_key(key):
                        self.AD_results[key] = {gene:{str(sorted_pos[i]):record}}
                        self.AR_dbHET_results[key] = {gene:{str(sorted_pos[i]):record}}
                    else:
                        if not self.AD_results[key].has_key(gene):
                            self.AD_results[key][gene] = {str(sorted_pos[i]):record}
                            self.AR_dbHET_results[key][gene] = {str(sorted_pos[i]):record}
                        else:
                            self.AD_results[key][gene][str(sorted_pos[i])] = record
                            self.AR_dbHET_results[key][gene][str(sorted_pos[i])] = record
                    #print key, "\t".join(record)

        
        # contintue for double HETs model
        self.compundHetCleaning()
        
    def compundHetCleaning(self):
        '''
        Remove genes with only one HET vairants        
        '''
        for key in self.AR_dbHET_results.keys():
            for gene in self.AR_dbHET_results[key].keys():
                if len(self.AR_dbHET_results[key][gene]) <= 1:
                    del self.AR_dbHET_results[key][gene]



def splitPaths(option, opt, value, parser):
    try:
        setattr(parser.values, option.dest, value.split(','))
    except Exception as e:
        print "For option -a you need at least paths to 2 VCF files sperated by space."
        print e
        sys.exit(0) 

def main():

    parser = optparse.OptionParser()
    parser.add_option('-a', '--affected',  type='string',help='path to affected single VCF file)', action="callback", callback=splitPaths)

    
    parser.add_option('-f', '--filter', help='path to filters.txt', default="filters.txt")
    parser.add_option('-r', '--hierarchy', help='path to hierarchy.txt', default="hierarchy.txt")
    parser.add_option('-o', '--output', help='path to output.txt', default="output.txt")
    parser.add_option('-t', '--tags', help='path to tags.txt', default="tags.txt")
    parser.add_option('-s', '--save', help='Excel or R', default=None)

    (opts, args) = parser.parse_args()

    #paths = ["../../EVAR/examples/sibs/sib1_5k.vcf.gz", "../../EVAR/examples/sibs/sib2_5k.vcf.gz","../../EVAR/examples/sibs/sib3_5k.vcf.gz"]
    #opts.affected = '/Users/Macia/Desktop/DropBox/Dropbox/Dropbox/Team29/MyProjects/CHD/Patel_Eammon/data/2011-11-11/Consanguineous/Merged_with_UK10Kmaf/D11_12780.vc_UK10Kmaf.vcf.gz /Users/Macia/Desktop/DropBox/Dropbox/Dropbox/Team29/MyProjects/CHD/Patel_Eammon/data/2011-11-11/Consanguineous/Merged_with_UK10Kmaf/D11_12781.vc_UK10Kmaf.vcf.gz /Users/Macia/Desktop/DropBox/Dropbox/Dropbox/Team29/MyProjects/CHD/Patel_Eammon/data/2011-11-11/Consanguineous/Merged_with_UK10Kmaf/D11_12784.vc_UK10Kmaf.vcf.gz'
    paths = opts.affected
    
    
    if None in paths or len(paths) < 1:
        print paths
        print "Please make sure you have at least 2 VCF files. Use option -h to print the help message"
        sys.exit(0)


    filters_path = opts.filter
    weights_path = opts.hierarchy
    output_path = opts.output
    tags_path  = opts.tags

    mySingle = Single(paths, weights_path, filters_path, output_path, tags_path, opts.save)

if __name__ == "__main__":
    main()
    

