# Write a strataG gtypes object from GDS or tidy data

Write a [strataG](https://github.com/EricArcher/strataG) object from a
tidy data frame. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
write_gtypes(data, write = FALSE, filename = NULL)
```

## Arguments

- data:

  A supported genomic file, object, or tidy genomic data frame. Default:
  `data = NULL`.

- write:

  (logical, optional) To write in the working directory the gtypes
  object. The file is written with `radiator_gtypes_DATE@TIME.RData` if
  no filename is provided and can be open with load or readRDS. Default:
  `write = FALSE`.

- filename:

  (character, optional) Filename prefix. Default: `filename = NULL`.

## Value

An object of the class [strataG](https://github.com/EricArcher/) is
returned.

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

`gtypes` support requires the optional GitHub package
[strataG](https://github.com/EricArcher/strataG). Consult its repository
for current installation and troubleshooting information.

## Data filtering

This writer does not silently filter markers or individuals. It may
validate requirements imposed by the destination format and stop with an
informative error when the input is unsuitable. It is the user's
responsibility to filter and quality-control the data appropriately for
the intended analysis before generating the output. Use
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow when filtering is required.

## References

Archer FI, Adams PE, Schneiders BB. strataG: An r package for
manipulating, summarizing and analysing population genetic data.
Molecular Ecology Resources. 2017; 17: 5-11. doi:10.1111/1755-0998.12559

## See also

[strataG](https://github.com/EricArcher/)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# require(strataG)
# with radiator GDS
turtle <- genometranslator::write_gtypes(data = "my.metadata.node.rad")

# with tidy data
turtle <- genometranslator::write_gtypes(data = "my.radiator.rad")
} # }
```
