# Read a genlight object into a tidy data frame and/or GDS object/file

Read a genlight object into a tidy data frame and/or GDS object/file.
Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
read_genlight(
  data,
  tidy = TRUE,
  gds = TRUE,
  write = FALSE,
  verbose = FALSE,
  parallel.core = parallel::detectCores() - 1
)
```

## Arguments

- data:

  (path or object) A genlight object in the global environment or path
  to a genlight file that will be open with `readRDS`.

- tidy:

  (logical) Generate a tidy dataset. Default: `tidy = TRUE`.

- gds:

  (optional, logical) To write a radiator gds object. Default:
  `gds = TRUE`.

- write:

  (optional, logical) To write in the working directory the tidy data.
  The file is written with `radiator_genlight_DATE@TIME.rad`. Default:
  `write = FALSE`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

## Details

A string of the same dimension is generated when genlight:

1.  `is.null(genlight@pop)`: pop will be integrated in the tidy dataset.

2.  `is.null(data@chromosome)`: CHROM1 will be integrated in the tidy
    dataset.

3.  `is.null(data@loc.names)`: LOCUS1 to `ncol(genlight)` will be
    integrated in the tidy dataset.

4.  `is.null(data@position)`: an integer string of length =
    `ncol(genlight)` will be integrated in the tidy dataset.

**Note: that if all CHROM, LOCUS and POS is missing the function will be
terminated**

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

Reading and writing `genlight` objects requires the optional CRAN
package [adegenet](https://cran.r-project.org/package=adegenet).

## References

Jombart T (2008) adegenet: a R package for the multivariate analysis of
genetic markers. Bioinformatics, 24, 1403-1405.

Jombart T, Ahmed I (2011) adegenet 1.3-1: new tools for the analysis of
genome-wide SNP data. Bioinformatics, 27, 3070-3071.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
