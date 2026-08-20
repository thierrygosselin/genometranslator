# Write a GDS object from a tidy data frame

Write a Genomic Data Structure (GDS) file format
[gdsfmt](https://github.com/zhengxwen/gdsfmt) and object of class
`SeqVarGDSClass` from a tidy data frame.

Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
write_gds(
  data,
  data.source = NULL,
  filename = NULL,
  open = TRUE,
  verbose = TRUE
)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- data.source:

  (optional, character) The name of the software that generated the
  data. e.g. `data.source = "Stacks v.2.4"`. Default:
  `data.source = NULL`.

- filename:

  (optional) The file name of the Genomic Data Structure (GDS) file.
  radiator will append `.gds.rad` to the filename. If filename chosen is
  already present in the working directory, the default
  `radiator_datetime.gds.rad` is chosen. Default: `filename = NULL`.

- open:

  (optional, logical) Open or not the radiator connection. Default:
  `open = TRUE`.

## Value

An object in the global environment of class `SeqVarGDSClass` and a file
in the working directory.

## Data filtering

This writer does not silently filter markers or individuals. It may
validate requirements imposed by the destination format and stop with an
informative error when the input is unsuitable. It is the user's
responsibility to filter and quality-control the data appropriately for
the intended analysis before generating the output. Use
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow when filtering is required.

## Dependencies

Required package dependencies are declared in `DESCRIPTION` and are
installed with genometranslator. Any additional dependency needed only
for this format or option is identified in this help page. Use
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect the availability of core packages, optional packages, and
external executables.

## References

Zheng X, Levine D, Shen J, Gogarten SM, Laurie C, Weir BS. A
high-performance computing toolset for relatedness and principal
component analysis of SNP data. Bioinformatics. 2012;28: 3326-3328.
doi:10.1093/bioinformatics/bts606

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C,
Levine D (2017). SeqArray – A storage-efficient high-performance data
format for WGS variant calls. Bioinformatics.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
data.gds <- genometranslator::write_gds(data = "shark.rad")
} # }
```
