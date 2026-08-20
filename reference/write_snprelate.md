# Write a SNPRelate object from a tidy data frame

Write a [SNPRelate](https://github.com/zhengxwen/SNPRelate) object from
a tidy data frame. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users. SNPRelate analyses commonly require
biallelic SNPs and may require linkage disequilibrium pruning. Make
those filtering decisions before export for the intended SNPRelate
analysis.

## Usage

``` r
write_snprelate(data, biallelic = TRUE, filename = NULL, verbose = TRUE)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).
  **The genotypes are biallelic.**

- biallelic:

  (logical, optional) If you already know that the data is biallelic use
  this argument to speed up the function. Default: `biallelic = TRUE`.

- filename:

  (optional) The file name of the Genomic Data Structure (GDS) file.
  radiator will append `.gds` to the filename. If filename chosen is
  already present in the working directory, the default
  `radiator_snprelate_datetime.gds` is chosen. Default:
  `filename = NULL`.

## Value

An object in the global environment of class
`"SNPGDSFileClass", "gds.class"` and a file in the working directory.

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

This writer requires the optional Bioconductor package
[SNPRelate](https://bioconductor.org/packages/SNPRelate).

## Data filtering

This writer does not silently filter markers or individuals. It may
validate requirements imposed by the destination format and stop with an
informative error when the input is unsuitable. It is the user's
responsibility to filter and quality-control the data appropriately for
the intended analysis before generating the output. Use
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow when filtering is required.

## References

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS. A
high-performance computing toolset for relatedness and principal
component analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
doi:10.1093/bioinformatics/bts606

## See also

[SNPRelate](https://github.com/zhengxwen/SNPRelate)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
require(SNPRelate)
data.gds <- genometranslator::write_snprelate(data = "shark.rad")
} # }
```
