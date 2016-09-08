[![Build Status](https://travis-ci.org/jeremymcrae/clinical-filter.svg?branch=master)]
(https://travis-ci.org/jeremymcrae/clinical-filter)
[![Coverage Status](https://coveralls.io/repos/github/jeremymcrae/clinical-filter/badge.svg?branch=master)]
(https://coveralls.io/github/jeremymcrae/clinical-filter?branch=master)

## Clinical filtering for trios
Find candidate diagnostic variants in affected children that might contribute to
their disorder. We load VCF files (either named on the command line, or listed
in a PED file) for members of a family, filter for rare, functionally disruptive
variants, and assess whether each variant might affect the child's disorder. We
take into account the parents genotypes (if available) and whether the parents
are also affected with a (the?) disorder. For variants in known disease
causative genes we check whether the inheritance patterns matches one expected
for the inheritance models of the gene.

### VCF requirements
The code expects VCFs as version 4.2. Many INFO entries need to take multiple
alleles, and/or multiple genes into account. Multi-allelic variants expect
comma-separated entries for many fields, such as in the HGNC field e.g
'GENE1,GENE1' for a multi-allelic variant where both alleles affect the gene
'GENE1'. Multi-genic variants similarly expect '|' separated entries for the
multiple genes e.g. 'GENE1|GENE2' for a variant that occurs in both 'GENE1' and
'GENE2'.

Consequence strings come from VEP, and are expected in the `CQ` entry in the
INFO. De novo mutations need a `PP_DNM` (posterior probability of de novo
mutation, estimated by denovogear) entry in the FORMAT. By default, we screen
for de novos with PP_DNM > 0.9.

### Install
```sh
pip install git+git://github.com/jeremymcrae/clinical-filter.git --user

# Alternatively:
git clone https://github.com/jeremymcrae/clinical-filter.git
cd clinical-filter
python setup.py install --user
```

### Usage
For running the filtering, the basic command is either with a ped file
specified, i.e.

```sh
python clinical_filter.py \
  --ped PED_PATH
```

Or the individual VCF files for a trio can be specified, i.e.

```sh
python clinical_filter.py \
  --child CHILD_VCF_PATH \
  --mother MOTHER_VCF_PATH \
  --father FATHER_VCF_PATH \
  --gender M \ #M|F
  --mom-aff MOM_AFFECTED_STATUS (1=unaffected or 2=affected) \
  --dad-aff DAD_AFFECTED_STATUS (1=unaffected or 2=affected)
```

The ped option is the easiest if you have a large number of trios to
process, so you can define all the families and their VCF paths in the ped
file, and run with that.

Other options are:
 * `--syndrome-regions SYNDROMES_PATH` # path to file listing DECIPHER regions
 * `--known-genes KNOWN_GENES_PATH` # to specify the DDG2P database file
 * `--known-genes-date 2014-01-01` # to specify the version of the known genes file
 * `--alternate-ids ALTERNATE_IDS_PATH` # path to file for mapping individuals
   between IDs used in the PED file, to alternate study IDs.
 * `--output OUTPUT_PATH` # to specify that you want tab-separated output
   written to the given path
 * `--export-vcf OUTPUT_VCF_PATH` # to specify you want a filtered VCF, can
   give a directory (when analysing multiple individuals), or give a file path
 * `--maf-populations POP1_AF,POP2_AF` # to specify populations with MAF values
   within the INFO field of variants.

The output options can be omitted, or used together, whichever you need.
