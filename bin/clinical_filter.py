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

def main():
    """ run the clinical filtering analyses
    """
    
    options = get_options()
    
    # set the level of logging to generate
    numeric_level = getattr(logging, options.loglevel.upper(), None)
    
    log_filename = "clinical-filter.log"
    if options.ped is not None:
        log_filename = options.ped + ".log"
    
    logging.basicConfig(level=numeric_level, filename=log_filename)
    
    finder = Filter(options)
    finder.filter_trios()

if __name__ == "__main__":
    main()
