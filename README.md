[![Build Status](https://travis-ci.org/jeremymcrae/clinical-filter.svg?branch=master)]
(https://travis-ci.org/jeremymcrae/clinical-filter)
[![Coverage Status](https://coveralls.io/repos/github/jeremymcrae/clinical-filter/badge.svg?branch=master)]
(https://coveralls.io/github/jeremymcrae/clinical-filter?branch=master)

## Clinical filtering for trios

Find variants in affected children that might contribute to their disorder. We
load VCF files (either named on the command line, or listed in a PED file) for
members of a family, filter for rare, functionally disruptive variants, and
assess whether each variant might affect the child's disorder. We take into
account the parents genotypes (if available) and whether the parents are also
affected with a (the?) disorder. For variants in known disease causative genes
we check whether the inheritance patterns matches one expected for the
inheritance models of the gene.

### VCF requirements
Gene symbols are expected as either `HGNC` (single gene) or `HGNC_ALL`
("&"-seperated multiple genes) entries in the INFO field. Consequence strings
come from VEP, and are expected in the `CQ` entry in the INFO. De novo mutations
need a `PP_DNM` (posterior probability of de novo mutation, estimated by
denovogear) entry in the FORMAT. By default, we screen for de novos with
PP_DNM > 0.9.

We only check rare variants by excluding variants where any reference population
has a minor allele frequency (MAF) >= 0.01. The MAF for the reference
populations are included in the INFO field, if the MAF is available for that
population. Currently the populations that are checked are:
* [continental 1000 Genomes populations](http://www.1000genomes.org/about)
* [DDD unaffected parents population](http://www.ddduk.org/)
* [UK10K population](http://www.uk10k.org/)

The reference population tags are hard-coded as a class variable in
`../clinicalfilter/variant/info.py`. The tags for these populations in the VCF
are:
- AFR_AF
- AMR_AF
- ASN_AF
- DDD_AF
- EAS_AF
- ESP_AF
- EUR_AF
- MAX_AF
- SAS_AF
- UK10K_cohort_AF

### Usage
```sh
python clinical_filter.py \
  --ped temp_name.ped \
  --syndrome-regions regions_filename.txt \ # optional
  --known-genes known_genes.txt \ # optional (but recommended)
  --known-genes-date 2014-01-01 \ # optional
  --alternate-ids alternate_ids.txt \ # optional
  --output output_name.txt \ # optional
  -export-vcf vcf_filename_or_directory # optional
```

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

The output options can be omitted, or used together, whichever you need.
