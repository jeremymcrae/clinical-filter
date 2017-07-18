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

import unittest

from clinicalfilter.variant.symbols import Symbols

class TestVariantSymbolsPy(unittest.TestCase):
    """  unit testing of the Symbols class
    """
    
    def setUp(self):
        info = {'HGNC': 'A|B,C|D', 'HGNC_ID': '1|2,3|', 'SYMBOL': 'Z|H,|'}
        self.symbols = Symbols(info, 0)
    
    def test___repr__(self):
        ''' test Symbols repr
        '''
        self.assertEqual(repr(self.symbols), "Symbols(info={'ENSG': '|', " \
            "'ENSP': '|', 'ENSR': '|', 'ENST': '|', 'HGNC': 'A|B', " \
            "'HGNC_ID': '1|2', 'SYMBOL': 'Z|H'}, idx=0)")
    
    def test_prioritise(self):
        ''' test that we correctly prioritise gene symbols
        '''
        self.assertEqual(self.symbols.prioritise(), ['1', '2'])
        self.assertEqual(self.symbols.prioritise(priority=['HGNC']), ['A', 'B'])
        self.assertEqual(self.symbols.prioritise(priority=['ENST', 'HGNC']), ['A', 'B'])
    
    def test_get_preferred(self):
        ''' tets that we can get a symbol, prioritising by symbol type
        '''
        values = {'HGNC': 'A', 'HGNC_ID': '1', 'SYMBOL': 'Z', 'ENSG': None,
            'ENST': None, 'ENSP': None, 'ENSR': None}
        
        # defaul to HGNC_ID first
        self.assertEqual(self.symbols.get_preferred(values), '1')
        
        # if we provide a list of symbols, check that order instead
        self.assertEqual(self.symbols.get_preferred(values, ['HGNC']), 'A')
        self.assertEqual(self.symbols.get_preferred(values, ['ENST']), None)
        
        # run through the list of preferred symbol types until we hit the end,
        # or get a non-None value
        self.assertEqual(self.symbols.get_preferred(values, ['ENST', 'HGNC']), 'A')
    
    def test_get(self):
        ''' test that we can retrieve gene symbols
        '''
        self.assertEqual(self.symbols.get('A'), '1')
        self.assertEqual(self.symbols.get('A', ['SYMBOL']), 'Z')
        
        self.assertEqual(self.symbols.get('A', 'ENST'), None)
        self.assertEqual(self.symbols.get('A', ['ENST']), None)
        self.assertEqual(self.symbols.get('A', ['ENST', 'SYMBOL']), 'Z')
        
        with self.assertRaises(KeyError):
            self.symbols.get('C', ['ENST'])
            self.symbols.get('A', 'UNKNOWN')
    
    def test_second_allele(self):
        ''' test that we can set gene symbols from the second allele
        '''
        info = {'HGNC': 'A|B,C|D', 'HGNC_ID': '1|2,3|', 'SYMBOL': 'Z|H,|'}
        symbols = Symbols(info, 1)
        
        self.assertEqual(symbols.get('C'), '3')
        
        # test that if a symbol is missing, we skip to an alternate
        self.assertEqual(symbols.get('D'), 'D')
        
        # check that we can't retrieve anything for a symbol in the other allele
        with self.assertRaises(KeyError):
            symbols.get('A')
    
    def test_out_of_index_allele(self):
        ''' raise an error if we construct a class for a non-exitent allele
        '''
        info = {'HGNC': 'A|B,C|D', 'HGNC_ID': '1|2,3|', 'SYMBOL': 'Z|H,|'}
        with self.assertRaises(IndexError):
            Symbols(info, 2)
    
    def test_set(self):
        ''' test we can set gene symbols after the class has been instantiated
        '''
        
        # try a key that does not currently exist
        self.assertEqual(self.symbols.get('A', 'ENST'), None)
        
        # set the key
        self.symbols.set('A', 'H', 'ENST')
        self.assertEqual(self.symbols.get('A', 'ENST'), 'H')
        
        # overwrite an existsing key
        self.symbols.set('A', 'O', 'HGNC')
        self.assertEqual(self.symbols.get('O', 'HGNC'), 'O')
        
        with self.assertRaises(KeyError):
            self.symbols.get('A')
        
