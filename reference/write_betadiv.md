# Write a betadiv file from a tidy data frame

Write a betadiv file from a tidy data frame. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and [assigner](https://github.com/thierrygosselin/assigner) and might be
of interest for users.

## Usage

``` r
write_betadiv(data)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

## Value

A betadiv object is returned.

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

Lamy T, Legendre P, Chancerelle Y, Siu G, Claudet J (2015) Understanding
the Spatio-Temporal Response of Coral Reef Fish Communities to Natural
Disturbances: Insights from Beta-Diversity Decomposition. PLoS ONE, 10,
e0138696.

## See also

`beta.div` is available on Pierre Legendre web site
<http://adn.biol.umontreal.ca/~numericalecology/Rcode/>

## Author

Laura Benestan <laura.benestan@icloud.com> and Thierry Gosselin
<thierrygosselin@icloud.com>
