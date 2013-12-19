"""Additional functions to load user specified values to search for variants that fit specified 
inheritance models.
"""

import os

# Secify the first and second pseudoautosomal regions on the X-chromosome, as well as the 
# X-transposed region on the X-chromosome
pseudoautosomal_regions = [(1,2699520), (154930290,155260560), (88456802,92375509)]

def getTrioModel(path):
    """Opens values that test whether variants fit specified inheritance models.

    Opens a tab-separated file that specifies whether variants in trios fit inheritance models. The 
    file has columns which cover 1) the inheritance model (eg dominant, recessive or recessive 
    compound heterozygous), 2) the child, mother and fathers VCF coded genotypes, coded like a 
    python tuple, eg (1,2,0), 3) whether the specified genotypes fit the inheritance model when the 
    variant is on an autosomal chromosome, 4) whether the specified genotypes fit the inheritance 
    model when the variant is on a X-chromosome, and the child is male, and 5) whether the 
    specified genotypes fit the inheritance model when the variant is on an X-chromosome, and the 
    child is female. 

    The last three columns are True/False values, so are eval'ed for conditional testing.

    Args:
        path: path to text file defining the trio inheritance model.

    Returns:
        A dictionary of genotypes for the trio under difference models if the genotypes fit the 
        model. For example, under complete penetrance the dictionary should be:

        {'dominant': {'autosomal': [], 'XChrMale': [], 'XChrFemale': [(1,1,0)]}, 
        recessive': {'autosomal': [(2,1,1)], 'XChrMale': [(2,1,0)], 'XChrFemale': []}, 
        'recessiveCompoundHet': {'autosomal': [((1,1,0),(1,0,1))], 'XChrMale': [],'XChrFemale': []}}

    Raises:
        IOError: An error when the inheritance model file path is not specified correctly.

    NOTE: This function might be more suited to being in user.py, since that is where all the 
    other functions that load user defined values reside.
    """

    if not os.path.exists(path):
        raise IOError("path to the trio inheritance model file does not exist.")

    chrGroups = ['autosomal','XChrMale','XChrFemale']
    models = ['dominant','recessive','recessiveCompoundHet']
    affected_status = ["parents_unaffected", "father_affected", "mother_affected", "both_parents_affected"]
    
    # prepare a dictionary with all the chrGroups by models having blank lists
    trio_inheritance_table = {}
    for group in chrGroups:
        for status in affected_status:
            for model in models:
                # make sure the dict contains the affected status and the chomosome model
                if model not in trio_inheritance_table:
                    trio_inheritance_table[model] = {}
                
                if status not in trio_inheritance_table[model]:
                    trio_inheritance_table[model][status] = {}
                
                trio_inheritance_table[model][status][group] = []
    
    # go through the inheritance model file, and when a genotype passes the inheritance model, add 
    # the genotypes that pass to the correct dictionary entry.
    f = open(path)
    for line in f:
        if not line.startswith("#"):
            line = line.strip().split('\t')
            model = line[0]
            affected_status = line[1]
            genotypeSet = eval(line[2])
            for i,includeGenotype in enumerate(line[3:5]):
                includeGenotype = eval(includeGenotype)
                if includeGenotype:
                    trio_inheritance_table[model][affected_status][chrGroups[i]].append(genotypeSet)
    return trio_inheritance_table


