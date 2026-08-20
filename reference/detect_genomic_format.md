# Used internally in radiator to detect the file format

Detect the file format of genomic data set.

## Usage

``` r
detect_genomic_format(data, guess = NULL)
```

## Arguments

- data:

  15 options for input: VCFs (SNPs or Haplotypes, to make the vcf
  population ready), plink (tped and bed), stacks haplotype file, genind
  (library(adegenet)), genlight (library(adegenet)), gtypes
  (library(strataG)), genepop, DArT, and a data frame in long/tidy or
  wide format. To verify that radiator detect your file format use
  `detect_genomic_format` (see example below). New addition to radiator:
  the Apache Parquet columnar storage file format that will replace the
  fst (Lightning Fast Serialiation) format. Documented in **Input
  genomic datasets** of
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- guess:

  (character) In development, guess faster the type of file. Default:
  `guess = NULL`.

## Value

One of these file format:

- tbl_df: for a data frame

- genind: for a genind object

- genlight: for a genlight object

- gtypes: for a gtypes object

- vcf.file: for a vcf file

- plink.tped.file: for a plink tped file

- plink.bed.file: for a plink bed file

- genepop.file: for a genepop file

- haplo.file: for a stacks haplotypes file

- fstat.file: for a fstat file

- dart: for a DArT file

- fst.file: for a file ending with .rad

- SeqVarGDSClass: for SeqArray GDS file.

- arrow parquet: for a Apache Parquet columnar storage file format, used
  by arrow R package.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
#To verify your file is detected by radiator as the correct format:
genometranslator::detect_genomic_format(data = "populations.snps.vcf")
} # }
```
