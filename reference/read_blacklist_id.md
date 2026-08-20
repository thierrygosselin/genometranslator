# read_blacklist_id

Read a file or object with blacklisted individuals. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
read_blacklist_id(blacklist.id = NULL, verbose = TRUE)
```

## Arguments

- blacklist.id:

  (optional, path or object) A blacklist file in the working directory
  or object in the global environment. The data frame as 1 column (named
  `INDIVIDUALS`) and is filled with the individual IDs The ids are
  cleaned with
  [`clean_ind_names`](https://thierrygosselin.github.io/genometranslator/reference/clean_ind_names.md)
  for separators, only `-` are tolerated. Duplicates are removed
  automatically. Default: `blacklist.id = NULL`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = TRUE`.

## Value

A tibble with column `INDIVIDUALS`.

## Dependencies

Required package dependencies are declared in `DESCRIPTION` and are
installed with genometranslator. Any additional dependency needed only
for this format or option is identified in this help page. Use
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect the availability of core packages, optional packages, and
external executables.

## Examples

``` r
if (FALSE) { # \dontrun{
bl <- genometranslator::read_blacklist_id("blacklist.tsv")
} # }
```
