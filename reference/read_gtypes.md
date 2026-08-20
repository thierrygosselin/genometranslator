# Read a gtypes object into a tidy data frame

Transform a [strataG gtypes](https://github.com/EricArcher/strataG)
object into a tidy data frame. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
read_gtypes(data, verbose = FALSE)
```

## Arguments

- data:

  A gtypes object (\>= v.2.0.2) in the global environment.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

`gtypes` support requires the optional GitHub package
[strataG](https://github.com/EricArcher/strataG). Consult its repository
for current installation and troubleshooting information.

## References

Archer FI, Adams PE, Schneiders BB. strataG: An r package for
manipulating, summarizing and analysing population genetic data.
Molecular Ecology Resources. 2017; 17: 5-11. doi:10.1111/1755-0998.12559

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
