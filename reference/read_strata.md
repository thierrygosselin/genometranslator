# read strata

Read a strata object or file. The strata file contains thes individual's
metadata, the stratification: e.g. the population id and/or the sampling
sites (see details). Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
read_strata(
  strata,
  pop.id = FALSE,
  pop.levels = NULL,
  pop.labels = NULL,
  pop.select = NULL,
  blacklist.id = NULL,
  keep.two = FALSE,
  path.folder = getwd(),
  filename = NULL,
  verbose = FALSE
)
```

## Arguments

- strata:

  (path or object) The strata file or object. Additional documentation
  is available in `read_strata`. Use that function to
  whitelist/blacklist populations/individuals. Option to set
  `pop.levels/pop.labels` is also available.

- pop.id:

  (logical) When `pop.id = TRUE`, the strata returns the stratification
  colname `POP_ID`. With the default, Returns `STRATA`. Default:
  `pop.id = FALSE`.

- pop.levels:

  (optional, string) This refers to the levels in a factor. In this
  case, the id of the pop. Use this argument to have the pop ordered
  your way instead of the default alphabetical or numerical order. e.g.
  `pop.levels = c("QUE", "ONT", "ALB")` instead of the default
  `pop.levels = c("ALB", "ONT", "QUE")`. White spaces in population
  names are replaced by underscore. Default: `pop.levels = NULL`.

- pop.labels:

  (optional, string) Use this argument to rename/relabel your pop or
  combine your pop. e.g. To combine `"QUE"` and `"ONT"` into a new pop
  called `"NEW"`: (1) First, define the levels for your pop with
  `pop.levels` argument: `pop.levels = c("QUE", "ONT", "ALB")`. (2)
  then, use `pop.labels` argument:
  `pop.labels = c("NEW", "NEW", "ALB")`. To rename `"QUE"` to `"TAS"`:
  `pop.labels = c("TAS", "ONT", "ALB")`. Default: `pop.labels = NULL`.
  White spaces in population names are replaced by underscore.

- pop.select:

  (optional, string) Selected list of populations for the analysis. e.g.
  `pop.select = c("QUE", "ONT")` to select `QUE` and `ONT` population
  samples (out of 20 pops). If `pop.labels` argument was used to rename
  the strata column, use the new names with `pop.select`. Default:
  `pop.select = NULL`.

- blacklist.id:

  (optional, path or object) A blacklist file in the working directory
  or object in the global environment. The data frame as 1 column (named
  `INDIVIDUALS`) and is filled with the individual IDs The ids are
  cleaned with
  [`clean_ind_names`](https://thierrygosselin.github.io/genometranslator/reference/clean_ind_names.md)
  for separators, only `-` are tolerated. Duplicates are removed
  automatically. Default: `blacklist.id = NULL`.

- keep.two:

  (optional, logical) The output is limited to 2 columns:
  `INDIVIDUALS, STRATA`. By default all the samples metadata is
  imported. Default: `keep.two = FALSE`.

- path.folder:

  (optional, path) If `!is.null(blacklist.id) || !is.null(pop.select)`,
  the modified strata is written by default in the working directory.
  Default: `path.folder = getwd()`.

- filename:

  (optional, character) If
  `!is.null(blacklist.id) || !is.null(pop.select)`, the modified strata
  is written by default in the working directory with date and time
  appended to `strata_radiator_filtered`, to make the file unique. If
  you plan on writing more than 1 strata file per minute, use this
  argument to supply the unique filename. When filename is not NULL, it
  will also trigger saving the strata to a file. Default:
  `filename = NULL`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.

## Value

**A list** with several components:

1.  \$strata

2.  \$pop.levels

3.  \$pop.labels

4.  \$pop.select

5.  \$blacklist.id

## Details

The strata file used in radiator is a tab delimited file with a minimum
of 2 columns headers (3 for DArT data users): `INDIVIDUALS` and
`STRATA`. If a `strata` file is specified with all file formats that
don't require it, the strata argument will have precedence on the
population groupings used internally in those file formats. For file
formats without population/strata groupings (e.g. vcf, haplotype files)
if no strata file is provided, 1 pop/strata grouping will automatically
be created. For vcf and haplotypes file, the strata can also be used as
a whitelist of id. Samples not in the strata file will be discarded from
the data set. The `STRATA` column can be any hierarchical grouping. To
create a strata file see
[`individuals2strata`](https://thierrygosselin.github.io/genometranslator/reference/individuals2strata.md).
If you have already run
[stacks](http://catchenlab.life.illinois.edu/stacks/) on your data, the
strata file is similar to a stacks *population map file*, make sure you
have the required column names (`INDIVIDUALS` and `STRATA`). The strata
column is cleaned of a white spaces that interfere with some packages or
codes: space is changed to an underscore `_`.

For DArT data see
[`read_dart`](https://thierrygosselin.github.io/genometranslator/reference/read_dart.md)

[example.strata.tsv](https://www.dropbox.com/s/g0vsek0dmtpxntt/example.strata.tsv?dl=0).

[example.dart.strata.tsv](https://www.dropbox.com/s/utq2h6o00v55kep/example.dart.strata.tsv?dl=0).

## VCF

VCF file users, not sure about the sample id inside your file ? See the
example in
[`extract_individuals_vcf`](https://thierrygosselin.github.io/genometranslator/reference/extract_individuals_vcf.md)

## DArT

DArT file users, not sure about the sample id inside your file ? See the
example in
[`extract_dart_target_id`](https://thierrygosselin.github.io/genometranslator/reference/extract_dart_target_id.md)

## Dependencies

Required package dependencies are declared in `DESCRIPTION` and are
installed with genometranslator. Any additional dependency needed only
for this format or option is identified in this help page. Use
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect the availability of core packages, optional packages, and
external executables.

## See also

[`summary_strata`](https://thierrygosselin.github.io/genometranslator/reference/summary_strata.md),
[`individuals2strata`](https://thierrygosselin.github.io/genometranslator/reference/individuals2strata.md),
[`change_pop_names`](https://thierrygosselin.github.io/genometranslator/reference/change_pop_names.md),
[`join_strata`](https://thierrygosselin.github.io/genometranslator/reference/join_strata.md),
[`generate_strata`](https://thierrygosselin.github.io/genometranslator/reference/generate_strata.md)

## Examples

``` r
if (FALSE) { # \dontrun{
strata.info <- genometranslator::read_strata(strata)

# the return object is a list with 5 objects:
names(strata.info)

# to get the strata
new.strata <- strata.info$strata

# if naything is changed from the original strata, a new strata file is
# generated automatically:

new.strata <- genometranslator::read_strata(
    strata = strata,
    blacklist.id = "blacklisted.ids.tsv"
    )

} # }
```
