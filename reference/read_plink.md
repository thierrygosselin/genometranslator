# Reads PLINK tped and bed files

The function reads PLINK tped and bed files. radiator prefers the use of
BED file. These files are converted to a connection SeqArray
[SeqArray](https://github.com/zhengxwen/SeqArray) GDS object/file of
class `SeqVarGDSClass` (Zheng et al. 2017). The Genomic Data Structure
(GDS) file format is detailed in
[gdsfmt](https://github.com/zhengxwen/gdsfmt).

Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
read_plink(
  data,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
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

- filename:

  (optional) The file name of the Genomic Data Structure (GDS) file.
  radiator will append `.gds.rad` to the filename. If the filename
  chosen exists in the working directory, the default
  `radiator_datetime.gds` is chosen. Default: `filename = NULL`.

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

- ...:

  Additional arguments passed to lower-level readers, translators, or
  writers.

## Value

For `tped` the function returns a list object with the non-modified
`tped` and the strata corresponding to the `tfam`. With `bed`, the
function returns a GDS object.

## Details

Large PLINK files will require the use of BED plink format. Look below
in the example for conversion with PLINK.

Large PLINK bed files will take longer to import and transform in GDS,
but after the file is generated, you can close your computer and come
back to it a month later and it's now a matter of sec to open a
connection!

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

Binary BED input uses the declared Bioconductor dependencies SeqArray
and gdsfmt. The external PLINK executable is not required to read or
write the supported PLINK formats.

## References

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C,
Levine D (2017). SeqArray – A storage-efficient high-performance data
format for WGS variant calls. Bioinformatics.

PLINK: a tool set for whole-genome association and population-based
linkage analyses. American Journal of Human Genetics. 2007: 81: 559–575.
doi:10.1086/519795

## See also

[PLINK](https://www.cog-genomics.org/plink/1.9/)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data <- genometranslator::read_plink(data = "my_plink_file.bed")
# when conversion is required from TPED to BED, in Terminal:
# plink --tfile my_plink_file --make-bed --allow-no-sex --allow-extra-chr --chr-set 95
} # }
```
