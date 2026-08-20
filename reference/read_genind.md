# Read a genind object to a GDS or tidy dataframe

Read a genind object or file from
[adegenet](https://github.com/thibautjombart/adegenet) to a tidy
dataframe. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
read_genind(data, tidy = TRUE, gds = TRUE, write = FALSE, verbose = FALSE)
```

## Arguments

- data:

  (path or object) A genind object in the global environment or path to
  a genind file that will be open with `readRDS`.

- tidy:

  (logical) Generate a tidy dataset. Default: `tidy = TRUE`.

- gds:

  (optional, logical) To write a radiator gds object. Currently, for
  biallelic datasets only. Default: `gds = TRUE`.

- write:

  (optional, logical) To write in the working directory the tidy data.
  The file is written with `radiator_genind_DATE@TIME.rad`. Default:
  `write = FALSE`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.

## Note

[genind](https://github.com/thibautjombart/adegenet) objects, like
genepop, are not optimal genomic format for RADseq datasets, they lack
important genotypes and markers metadata: chromosome, locus, snp,
position, read depth, allele depth, etc.
[genlight](https://github.com/thibautjombart/adegenet) object is a more
interesting container and is memory efficient, see
[`read_genlight`](https://thierrygosselin.github.io/genometranslator/reference/read_genlight.md).

By default allele names will be kept for the tidy dataset, if the
alleles is numeric and length \< 3.

In the unlikely event that the genind object as no
stratification/population, *pop* will be added to the strata column.

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

Reading and writing `genind` objects requires the optional CRAN package
[adegenet](https://cran.r-project.org/package=adegenet).

## References

Jombart T (2008) adegenet: a R package for the multivariate analysis of
genetic markers. Bioinformatics, 24, 1403-1405.

Jombart T, Ahmed I (2011) adegenet 1.3-1: new tools for the analysis of
genome-wide SNP data. Bioinformatics, 27, 3070-3071.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
