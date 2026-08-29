# Extract individuals from vcf file

Function that returns the individuals present in a vcf file. Useful to
create a strata file or to make sure you have the right individuals in
your VCF. Only the VCF header and first variant record are read;
variants are not counted or imported. The fixed VCF columns and record
width are validated before sample IDs are returned.

## Usage

``` r
extract_individuals_vcf(vcf)
```

## Arguments

- vcf:

  (character, path) The path to the vcf file.

## Value

A tibble with a column: `INDIVIDUALS`.

## See also

genometranslator
[`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Built a strata file:
strata <- genometranslator::extract_individuals_vcf("my.vcf") %>%
    dplyr::mutate(STRATA = "fill this") %>%
    readr::write_tsv(x = ., file = "my.new.vcf.strata.tsv")
} # }
```
