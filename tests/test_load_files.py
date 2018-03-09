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

import tempfile
import shutil

import unittest
from clinicalfilter.load_files import get_header_positions, \
    parse_gene_line,  open_known_genes, open_cnv_regions, open_x_lr2_file

class TestLoadFilesPy(unittest.TestCase):
    ''' test the file loading functions
    '''
    
    @classmethod
    def setUpClass(cls):
        cls.tempdir = tempfile.mkdtemp()
    
    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tempdir)
    
    def setUp(self):
        ''' make sure we are only going to be working on temporary files
        '''
        
        self.temp = tempfile.NamedTemporaryFile(dir=self.tempdir, delete=False)
    
    def test_get_header_positions(self):
        ''' check that get_header_positions() works correctly
        '''
        
        self.temp.write('chrom\tstart\tstop\thgnc\n'.encode('utf8'))
        self.temp.flush()
        
        with open(self.temp.name) as handle:
            self.assertEqual(get_header_positions(handle, ['chrom', 'stop',
                'hgnc']), {'chrom': 0, 'stop': 2, 'hgnc': 3})
        
        # raise an error if we try to get a column name that does not appear in
        # the header
        with open(self.temp.name) as handle:
            with self.assertRaises(ValueError):
                get_header_positions(handle, ['chrom', 'stop', 'missing'])
    
    def test_parse_gene_line(self):
        ''' check that parse_gene_line() works correctly
        '''
        
        header = {'gene': 0, 'chr': 1, 'start': 2, 'stop': 3, 'type': 4,
            'mode': 5, 'mech': 6, 'hgnc_id': 7}
        
        line = ['TEST', 'chr1', '1000', '2000', 'confirmed dd gene',
            'Biallelic', 'Loss-of-function', '1001']
        
        self.assertEqual(parse_gene_line(line, header), ('1001', {
            'chrom': 'chr1', 'start': 1000, 'end': 2000,
            'symbol': 'TEST', 'status': set(['confirmed dd gene']),
            'inh': {'Biallelic': set(['Loss-of-function'])}}))
    
    def test_parse_gene_line_both_mechanism(self):
        ''' check that parse_gene_line() works correctly for 'Both' mechanism.
        '''
        
        header = {'gene': 0, 'chr': 1, 'start': 2, 'stop': 3, 'type': 4,
            'mode': 5, 'mech': 6, 'hgnc_id': 7}
        line = ['TEST', 'chr1', '1000', '2000', 'confirmed dd gene',
            'Both', 'Loss-of-function', '1001']
        
        self.assertEqual(parse_gene_line(line, header), ('1001', {
            'chrom': 'chr1', 'start': 1000, 'end': 2000,
            'symbol': 'TEST', 'status': set(['confirmed dd gene']),
            'inh': {'Biallelic': set(['Loss-of-function']),
                'Monoallelic': set(['Loss-of-function']),
                'Both': set(['Loss-of-function'])}}))
    
    def test_open_known_genes(self):
        ''' test that open_known_genes() works correctly
        '''
        
        header = ['gene', 'chr', 'start', 'stop', 'type', 'mode', 'mech', 'hgnc_id']
        line = ['TEST', '1', '1000', '2000', 'confirmed dd gene',
            'Monoallelic', 'Loss-of-function', '1001']
        
        self.temp.write(('\t'.join(header) + '\n').encode('utf8'))
        self.temp.write(('\t'.join(line) + '\n').encode('utf8'))
        self.temp.flush()
        
        self.assertEqual(open_known_genes(self.temp.name),
            {'1001': {'chrom': '1', 'start': 1000, 'end': 2000,
                'symbol': 'TEST', 'status': set(['confirmed dd gene']),
                'inh': {'Monoallelic': set(['Loss-of-function'])}}
            })
    
    def test_open_known_genes_multimodes(self):
        ''' test that open_known_genes() works correctly for genes with >1 modes
        '''
        
        header = ['gene', 'chr', 'start', 'stop', 'type', 'mode', 'mech', 'hgnc_id']
        line1 = ['TEST', '1', '1000', '2000', 'confirmed dd gene',
            'Monoallelic', 'Loss-of-function', '1001']
        line2 = ['TEST', '1', '1000', '2000', 'confirmed dd gene',
            'Biallelic', 'Loss-of-function', '1001']
        
        self.temp.write(('\t'.join(header) + '\n').encode('utf8'))
        self.temp.write(('\t'.join(line1) + '\n').encode('utf8'))
        self.temp.write(('\t'.join(line2) + '\n').encode('utf8'))
        self.temp.flush()
        
        self.assertEqual(open_known_genes(self.temp.name),
            {'1001': {'chrom': '1', 'start': 1000, 'end': 2000,
                'symbol': 'TEST', 'status': set(['confirmed dd gene']),
                'inh': {'Monoallelic': set(['Loss-of-function']),
                     'Biallelic': set(['Loss-of-function'])}}
            })
    
    def test_open_known_genes_multimechs(self):
        ''' test that open_known_genes() works correctly for genes with >1 mechs
        '''
        
        header = ['gene', 'chr', 'start', 'stop', 'type', 'mode', 'mech', 'hgnc_id']
        line1 = ['TEST', '1', '1000', '2000', 'confirmed dd gene',
            'Monoallelic', 'Loss-of-function', '1001']
        line2 = ['TEST', '1', '1000', '2000', 'confirmed dd gene',
            'Monoallelic', 'Activating', '1001']
        
        self.temp.write(('\t'.join(header) + '\n').encode('utf8'))
        self.temp.write(('\t'.join(line1) + '\n').encode('utf8'))
        self.temp.write(('\t'.join(line2) + '\n').encode('utf8'))
        self.temp.flush()
        
        self.assertEqual(open_known_genes(self.temp.name),
            {'1001': {'chrom': '1', 'start': 1000, 'end': 2000,
                'symbol': 'TEST', 'status': set(['confirmed dd gene']),
                'inh': {'Monoallelic': set(['Loss-of-function', 'Activating'])}}
            })
    
    def test_open_known_genes_multigenes(self):
        ''' test that open_known_genes() works correctly for multiple genes
        '''
        
        header = ['gene', 'chr', 'start', 'stop', 'type', 'mode', 'mech', 'hgnc_id']
        line1 = ['TEST', '1', '1000', '2000', 'confirmed dd gene',
            'Monoallelic', 'Loss-of-function', '1001']
        line2 = ['TEST2', '1', '3000', '4000', 'confirmed dd gene',
            'Monoallelic', 'Loss-of-function', '2001']
        
        self.temp.write(('\t'.join(header) + '\n').encode('utf8'))
        self.temp.write(('\t'.join(line1) + '\n').encode('utf8'))
        self.temp.write(('\t'.join(line2) + '\n').encode('utf8'))
        self.temp.flush()
        
        self.assertEqual(open_known_genes(self.temp.name),
            {'1001': {'chrom': '1', 'start': 1000, 'end': 2000,
                'symbol': 'TEST', 'status': set(['confirmed dd gene']),
                'inh': {'Monoallelic': set(['Loss-of-function'])}},
            '2001': {'chrom': '1', 'start': 3000, 'end': 4000,
                'symbol': 'TEST2', 'status': set(['confirmed dd gene']),
                'inh': {'Monoallelic': set(['Loss-of-function'])}}
            })
    
    def test_open_known_genes_wrong_status(self):
        ''' test that open_known_genes() filters out genes without a good status
        '''
        
        header = ['gene', 'chr', 'start', 'stop', 'type', 'mode', 'mech', 'hgnc_id']
        line1 = ['TEST', '1', '1000', '2000', 'possible dd gene',
            'Monoallelic', 'Loss-of-function', '1001']
        line2 = ['TEST2', '1', '3000', '4000', 'confirmed dd gene',
            'Monoallelic', 'Loss-of-function', '2001']
        
        self.temp.write(('\t'.join(header) + '\n').encode('utf8'))
        self.temp.write(('\t'.join(line1) + '\n').encode('utf8'))
        self.temp.write(('\t'.join(line2) + '\n').encode('utf8'))
        self.temp.flush()
        
        self.assertEqual(open_known_genes(self.temp.name),
            {'2001': {'chrom': '1', 'start': 3000, 'end': 4000,
                'symbol': 'TEST2', 'status': set(['confirmed dd gene']),
                'inh': {'Monoallelic': set(['Loss-of-function'])}}
            })
    
    def test_open_known_genes_missing_lines(self):
        ''' test that open_known_genes() works correctly when we can't find any genes
        '''
        
        header = ['gene', 'chr', 'start', 'stop', 'type', 'mode', 'mech', 'hgnc_id']
        
        self.temp.write(('\t'.join(header) + '\n').encode('utf8'))
        self.temp.flush()
        
        # if we have checked the file, and there aren't any genes in it, this
        # raises an error, since the most likely explanation is that something
        # has gone wrong with the data file, and likely the line-endings
        with self.assertRaises(ValueError):
            open_known_genes(self.temp.name)
    
    def test_open_cnv_regions(self):
        ''' test that open_cnv_regions() works correctly
        '''
        
        lines = ['id_syndrome_feature\tid_syndrome\tcopy_number\tchr_start\tchr_end\tchr\n',
            '20\t1\t1\t1569197\t2110236\t4\tNA\t2650330\t149066\t1t\n']
        lines = [ x.encode('utf8') for x in lines ]
        
        self.temp.writelines(lines)
        self.temp.flush()
        
        self.assertEqual(open_cnv_regions(self.temp.name),
            {('4', '1569197', '2110236'): '1'})

    def test_open_x_lr2_file(self):
        '''test that open_x_lr2_file() works correct
        '''

        header = ['Proband', 'sum_X_l2r']
        line1 = ['DDDP012345', '1000']
        self.temp.write(('\t'.join(header) + '\n').encode('utf8'))
        self.temp.write(('\t'.join(line1) + '\n').encode('utf8'))
        self.temp.flush()

        self.assertEqual(open_x_lr2_file(self.temp.name),{'DDDP012345': '1000'})
    
