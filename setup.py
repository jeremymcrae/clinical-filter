"""
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
"""

from setuptools import setup

setup(
    name = "clinical-filter",
    version = "0.4.8",
    author = "Jeremy McRae",
    author_email = "jeremy.mcrae@sanger.ac.uk",
    description="Clinical filtering for trios.",
    long_description = ("Find variants in affected children that might  "
        "contribute to their disorder. We load VCF files (either named on the "
        "command line, or listed in a PED file) for members of a family, filter "
        "for rare, functionally disruptive variants, and assess whether each "
        "variant might affect the childs disorder. We take into account the "
        "parents genotypes (if available) and whether the parents are also "
        "affected with a (the?) disorder. For variants in known disease "
        "causative genes we check whether the inheritance patterns matches one "
        "expected for the inheritance models of the gene.\n"
        "Written by Jeremy McRae (jm33@sanger.ac.uk), derived from code by "
        "Saeed Al Turki and Jeff Barrett."),
    license = "MIT",
    packages=["clinicalfilter", 'clinicalfilter.variant'],
    install_requires=['pytabix >= 0.0.2',
    ],
    url='https://github.com/jeremymcrae/clinical-filter',
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
    ],
    test_suite="tests"
)
