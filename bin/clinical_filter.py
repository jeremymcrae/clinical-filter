'''
Copyright (c) 2016 Genome Research Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
'''

import logging

from clinicalfilter.load_options import get_options
from clinicalfilter.filter import Filter
from clinicalfilter.ped import load_families, Family

def get_families(args):
    """ loads a list of Family objects for multiple families, or a single trio
    """
    
    if args.ped is None:
        fam_id = 'blank_family_ID'
        family = Family(fam_id)
        family.add_child('child', args.mother, args.father, args.gender, '2', args.child)
        if args.mother is not None:
            family.add_mother('mother', '0', '0', '2',  args.mom_aff, args.mother)
        if args.father is not None:
            family.add_father('father',  '0', '0', '1', args.dad_aff, args.father)
        
        families = [family]
    else:
        families = load_families(args.ped)
    
    return families

def main():
    """ run the clinical filtering analyses
    """
    
    args = get_options()
    
    # set the level of logging to generate
    numeric_level = getattr(logging, args.loglevel.upper(), None)
    
    log_filename = "clinical-filter.log"
    if args.ped is not None:
        log_filename = args.ped + ".log"
    
    logging.basicConfig(level=numeric_level, filename=log_filename)
    
    families = get_families(args)
    count = sum([ y.is_affected() for x in families for y in x.children ])
    
    finder = Filter(args.populations, count, args.known_genes, args.genes_date,
        args.alternate_ids, args.regions, args.lof_sites, args.pp_filter,
        args.output, args.export_vcf, args.debug_chrom, args.debug_pos)
    
    for family in families:
        finder.filter_trio(family)

if __name__ == "__main__":
    main()
