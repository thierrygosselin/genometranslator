# Write a LDna object from a tidy data frame

Write a [LDna](https://github.com/petrikemppainen/LDna) object from a
biallelic tidy data frame. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
write_ldna(
  data,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  ...
)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).
  **The genotypes are biallelic.**

- filename:

  (optional) The file name of the LDna lower matrix file. Radiator will
  append `.ldna.rds` to the filename. If filename chosen is already
  present in the working directory, the default
  `radiator_datetime.ldna.rds` is chosen. With default,
  `filename = NULL`, no file is generated, only an object in the Global
  Environment. To read the data back into R, use
  readRDS("filename.ldna.rds"). Default: `filename = NULL`.

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.

- ...:

  (optional) To pass further argument for fine-tuning the function (see
  details).

## Details

The function requires
[SNPRelate](https://github.com/zhengxwen/SNPRelate) to prepare the data
for LDna.

To install SNPRelate: install.packages("BiocManager")
BiocManager::install("SNPRelate")

To install LDna: devtools::install_github("petrikemppainen/LDna")

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

This writer uses the optional Bioconductor package
[SNPRelate](https://bioconductor.org/packages/SNPRelate) to calculate
the LD matrix. It writes data for LDna; it does not install or run the
separate LDna software.

## Data filtering

This writer does not silently filter markers or individuals. It may
validate requirements imposed by the destination format and stop with an
informative error when the input is unsuitable. It is the user's
responsibility to filter and quality-control the data appropriately for
the intended analysis before generating the output. Use
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow when filtering is required.

## References

Kemppainen P, Knight CG, Sarma DK et al. (2015) Linkage disequilibrium
network analysis (LDna) gives a global view of chromosomal inversions,
local adaptation and geographic structure. Molecular Ecology Resources,
15, 1031-1045.

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS. (2012) A
high-performance computing toolset for relatedness and principal
component analysis of SNP data. Bioinformatics. 28: 3326-3328.
doi:10.1093/bioinformatics/bts606

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
