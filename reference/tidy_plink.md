# Tidy PLINK tped and bed files

Transform bi-allelic PLINK files in `.tped` or `.bed` formats into a
tidy dataset.

Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and [assigner](https://github.com/thierrygosselin/assigner) and might be
of interest for users. Ensure marker type, missingness, allele coding,
and any linkage pruning are appropriate for the intended PLINK analysis
before export.

## Usage

``` r
tidy_plink(
  data,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  The PLINK file.

  - bi-allelic data only. For haplotypes use VCF.

  - `tped` file format: the corresponding `tfam` file must be in the
    directory.

  - `bed` file format: IS THE PREFERRED format, the corresponding `fam`
    and `bim` files must be in the directory.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.

- ...:

  Additional arguments passed to lower-level readers, translators, or
  writers.

## Value

A tidy tibble of the PLINK file.

## Advance mode

*dots-dots-dots ...* allows to pass several arguments for fine-tuning
the function:

1.  `calibrate.alleles`: logical. For `tped` files, if
    `calibrate.alleles = FALSE` the function runs faster but REF/ALT
    alleles may not be calibrated. The default assumes the users or
    sotware producing the PLINK file calibrated the alleles. Default:
    `calibrate.alleles = FALSE`.

## References

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C,
Levine D (2017). SeqArray – A storage-efficient high-performance data
format for WGS variant calls. Bioinformatics.

PLINK: a tool set for whole-genome association and population-based
linkage analyses. American Journal of Human Genetics. 2007: 81: 559–575.
doi:10.1086/519795

## See also

[PLINK](https://www.cog-genomics.org/plink/1.9/)

[`read_plink`](https://thierrygosselin.github.io/genometranslator/reference/read_plink.md)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data <- genometranslator::tidy_plink(data = "my_plink_file.bed", verbose = TRUE)


# when conversion is required from TPED to BED, in Terminal:
# plink --tfile my_plink_file --make-bed --allow-no-sex --allow-extra-chr --chr-set 95
} # }
```
