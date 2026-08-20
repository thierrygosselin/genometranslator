# Write a fineRADstructure file from a tidy data frame

Write a [fineRADstructure](https://github.com/millanek/fineRADstructure)
file from a tidy data frame. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
write_fineradstructure(data, strata = NULL, filename = NULL)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- strata:

  (path or object) The strata file or object. Additional documentation
  is available in
  [`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md).
  Use that function to whitelist/blacklist populations/individuals.
  Option to set `pop.levels/pop.labels` is also available. Default:
  `strata = NULL`.

- filename:

  (optional) The file name prefix for the fineRADstructure file written
  to the working directory. With default: `filename = NULL`, the date
  and time is appended to `radiator_fineradstructure_`. Default:
  `filename = NULL`.

## Value

A fineRADstructure file is written in the working directory. An object
is also returned in the global environment.

## Details

fineRADstructure requires special formating for populations, it needs to
be ONLY LETTERS, not numbers. This information is then merged with the
sample id (converted to integers). radiator will generate a dictionary
file.

## Life cycle

It become increasingly difficult for me to follow all the different
naming schemes researcher uses, if they're is any strategy...
Consequently, I have abandoned the idea of formating with letters the
populations to generate the fineRADstrucure file. If you get an error
see the the details section.

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

Malinsky, M., Trucchi, E., Lawson, D., Falush, D. (2018). RADpainter and
fineRADstructure: population inference from RADseq data. Mol. Biol.
Evol. 35(5), 1284-1290. https://dx.doi.org/10.1093/molbev/msy023

## See also

[fineRADstructure](https://github.com/millanek/fineRADstructure)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
