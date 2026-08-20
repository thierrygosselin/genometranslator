# Calibrate REF and ALT alleles based on count

Calibrate REF and ALT alleles based on counts. The REF allele is
designated as the allele with more counts in the dataset. The function
will generate a REF and ALT columns.

**reference genome**: for people using a reference genome, the reference
allele terminology is different and is not based on counts...

Used internally in genometranslator and might be of interest for users.

## Usage

``` r
calibrate_alleles(
  data,
  biallelic = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  A genomic data set in the global environment tidy formats. See details
  for more info.

- biallelic:

  (optional) If `biallelic = TRUE/FALSE` will be use during multiallelic
  REF/ALT decision and speed up computations. Default:
  `biallelic = NULL`.

- parallel.core:

  (optional) The number of core used for parallel execution. This is no
  longer used. The code is as fast as it can. Using more cores will
  reduce the speed. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  (optional, logical) `verbose = TRUE` to be chatty during execution.
  Default: `verbose = FALSE`.

- ...:

  (optional) To pass further argument for fine-tuning the tidying
  (details below).

## Value

Depending if the input file is biallelic or multiallelic, the function
will output additional to REF and ALT column several genotype codings:

- `GT`: the genotype in 6 digits format with integers.

- `GT_VCF`: the genotype in VCF format with integers.

- `GT_VCF_NUC`: the genotype in VCF format with letters corresponding to
  nucleotide.

- `ALT_DOSAGE`: alternate-allele dosage for biallelic markers: `0`, `1`,
  or `2` copies of the ALT allele, and `NA` for missing genotypes.

## Details

**Input data:** A minimum of 4 columns are required (the rest are
considered metata info):

1.  `MARKERS`

2.  `POP_ID`

3.  `INDIVIDUALS`

4.  `GT` and/or `GT_VCF_NUC` and/or `GT_VCF`

*How to get a tidy data frame ?* genometranslator
[`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
