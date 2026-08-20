# Write a genepopedit flatten object

Write a genepopedit object from a tidy data frame or GDS file/object.

Why not use `genepopedit::genepop_flatten`?

- genepopedit requires a [specific type of
  genepop](https://github.com/rystanley/genepopedit#genepopedit), so if
  you don't want to manipulate your genepop file, radiator is an
  alternative.

- radiator follows guidelines highlighted here: [genepop
  format](http://genepop.curtin.edu.au/help_input.md), but the 3
  functions in radiator that reads genepop files:

  1.  [genome_translator](https://thierrygosselin.github.io/radiator/reference/genome_translator.html)

  2.  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md)

  3.  the underlying module:
      [read_genepop](https://thierrygosselin.github.io/genometranslator/reference/read_genepop.html)

  imports a *larger variety of genepop alternatives*, similarly to
  [adegenet](https://github.com/thibautjombart/adegenet) `read.genepop`
  function, only faster.

## Usage

``` r
write_genepopedit(data)
```

## Arguments

- data:

  A supported genomic file, object, or tidy genomic data frame. Default:
  `data = NULL`.

## Value

A genepopedit object in the global environment.

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

Stanley RRE, Jeffery NW, Wringe BF, DiBacco C, Bradbury IR (2017)
genepopedit: a simple and flexible tool for manipulating multilocus
molecular data in R. Molecular Ecology Resources, 17, 12-18.

## See also

[genepopedit](https://github.com/rystanley/genepopedit)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
