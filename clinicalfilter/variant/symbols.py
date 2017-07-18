'''
Copyright (c) 2017 Genome Research Ltd.

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

class Symbols(object):
    ''' represent gene symbols for an alt allele
    '''
    
    # symbol types, in the default preferred order
    fields = ["HGNC_ID", "HGNC", "SYMBOL", "ENSG", "ENST", "ENSP", "ENSR"]
    
    def __init__(self, info, idx):
        ''' initialise the object with all the symbols for an alt allele
        
        Args:
            info: info dictionary for a VCF variant. Must contain entries for
                HGNC_ID, HGNC, SYMBOL, ENSG, ENST, ENSP, and ENSR. These entries
                are comma-separated lists of symbol lists for the alt alleles.
                The list is ordered as per the alt alleles of the variant. Each
                allele entry is a pipe-separated list of symbols e.g. 'A|B,A|B'
            idx: index position for an alt allele
        '''
        
        # get lists of symbols for each symbol type, for one alt allele
        temp = []
        for x in self.fields:
            try:
                temp.append(info[x].split(",")[idx].split("|"))
            except KeyError:
                # if the field is not present, just include an empty list
                temp.append([])
        
        # replace missing symbol values with None
        fix_missing = lambda x: [ y if y not in ['.', ''] else None for y in x ]
        temp = [ fix_missing(x) for x in temp ]
        
        # make sure all of the symbol lists have the same length. Occasionally
        # we get a variant with differing lengths. The one example I've seen had
        # ENSR=,.|ENSR00000215586. This is first split by ',' (giving ['',
        # '.|ENSR00000215586']), then an entry is selected  and split by '|'.
        # Only the seond entry contains '|', so the lengths are discrepant.
        k = max(( len(x) for x in temp ))
        temp = [ x if len(x) == k else [None] * k for x in temp ]
        
        # swap the data to a list of dictionaries
        self.symbols = []
        for i in range(k):
            data = dict(zip(self.fields, ( temp[x][i] for x in range(len(self.fields)) )))
            self.symbols.append(data)
    
    def __repr__(self):
        info = {}
        for field in self.fields:
            values = ( x[field] for x in self.symbols )
            values = [ x if x is not None else '' for x in values ]
            info[field] = '|'.join(values)
        
        return 'Symbols(info={}, idx={})'.format(info, 0)
    
    def __eq__(self, other):
        
        if len(self.symbols) != len(other.symbols):
            return False
        
        return all(( self.symbols[x] == other.symbols[x] for x in range(len(self.symbols)) ))
    
    def prioritise(self, priority=None):
        ''' return gene symbols, giving priority to HGNC IDs vs ENST symbols
        
        Args:
            priority: list of symbol types e.g. ['HGNC', 'ENSG'], in a given
                priority order.
        '''
        
        return [ self.get_preferred(x, priority) for x in self.symbols ]
    
    def get_preferred(self, symbols, priority=None):
        ''' return a symbol, prioritising by symbol type
        
        Args:
            symbols: dictionary of possible symbols for a gene
            priority: list of symbol types e.g. ['HGNC', 'ENSG'], in a given
                priority order.
        '''
        
        if priority is None:
            priority = self.fields
        
        for field in priority:
            value = symbols[field]
            
            if value is not None:
                break
        
        return value
    
    def get(self, symbol, priority=None):
        ''' get a symbol, given a different alternate symbol
        
        I want to be able to look up the gene that contains a given symbol, and
        pick out a preferred alternate symbol.
        
        Args:
            symbol: required symbol to match on. Could be a HGNC symbol, HGNC ID,
                or ENSG symbol etc.
            priority: preferred symbol type e.g. 'HGNC', or list of types e.g.
                ['HGNC', 'SYMBOL', 'ENSG']. The code only checks these types,
                unless priority=None, then it checks all types.
        '''
        
        if priority is None:
            priority = self.fields
        
        if type(priority) == str:
            priority = [priority]
        
        for x in self.symbols:
            if symbol not in x.values():
                continue
            
            return self.get_preferred(x, priority)
        
        raise KeyError('{} not found in symbols'.format(symbol))
    
    def set(self, symbol, alternate, field):
        ''' update a symbol.
        
        Currently, this is only used to remove HGNC_ID symbols from CNVs which
        appear to affect a known DD gene, but which does not actually overlap
        the gene range.
        '''
        
        for x in self.symbols:
            if symbol not in x.values():
                continue
            
            x[field] = alternate
