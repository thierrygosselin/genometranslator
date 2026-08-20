# Write a rubias object

Write a rubias object from a tidy data frame or GDS file/object.

## Usage

``` r
write_rubias(
  data,
  strata = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
)
```

## Arguments

- data:

  A supported genomic file, object, or tidy genomic data frame. Default:
  `data = NULL`.

- strata:

  (optional, tibble file or object) This tibble of individual's metadata
  must contain four columns:
  `SAMPLE_TYPE, REPUNIT, COLLECTION, INDIVIDUALS`. Those columns are
  described in **rubias**. With the default, With default, `SAMPLE_TYPE`
  is filled with `reference`. `REPUNIT` and `COLLECTION` will be filled
  by the `STRATA` or `POP_ID` column found in the data. Default:
  `strata = NULL`.

- filename:

  The prefix for the name of the file written to the working directory.
  Default: `filename = NULL`. With default, only the rubias object is
  generated. The filename will be appended `_rubias.tsv`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

## Value

A rubias object in the global environment and a file is written in the
working directory if `filename` argument was used.

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

Moran BM, and Anderson, E. C. 2018 Bayesian inference from the
conditional genetic stock identification model. Can. J. Fish. Aquat.
Sci. 76: 551-560.

Anderson, Eric C., Robin S. Waples, and Steven T. Kalinowski. (2008) An
improved method for predicting the accuracy of genetic stock
identification. Canadian Journal of Fisheries and Aquatic Sciences 65,
7:1475-1486.

Anderson, E. C. (2010) Assessing the power of informative subsets of
loci for population assignment: standard methods are upwardly biased.
Molecular ecology resources 10, 4:701-710.

## See also

[rubias](https://github.com/eriqande/rubias): genetic stock
identification (GSI) in the tidyverse.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
