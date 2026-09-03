
# genometranslator <a href="https://thierrygosselin.github.io/genometranslator/"><img src="man/figures/logo.png" align="right" height="160" /></a>

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![packageversion](https://img.shields.io/badge/Package%20version-0.0.0.9000-orange.svg)](commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2026--09--03-brightgreen.svg)](/commits/master)
<!-- badges: end -->

`genometranslator` reads, standardizes, and writes individual genomic
data. The easiest interface is `read_genome()` and `write_genome()`; use
a format-specific `read_*()` or `write_*()` function when more control
is needed. Large genomic datasets remain in a file-backed GDS rather
than being expanded into an R data frame unless that conversion is
requested explicitly.

## Quick start

`read_genome()` detects VCF, DArT, PLINK, FSTAT, Genepop, GDS, and
supported R genomic objects and dispatches to the corresponding reader:

``` r
genome <- genometranslator::read_genome(
  data = "individuals.vcf.gz",
  strata = "strata.tsv"
)
```

For VCF-specific controls, call `read_vcf()` directly. Both interfaces
return an open SeqArray GDS object:

``` r
genome <- genometranslator::read_vcf(
  data = "individuals.vcf.gz",
  strata = "strata.tsv",
  filename = "individuals"
)
```

Use `tidy_genome()` only when an in-memory table is genuinely needed:

``` r
genotypes <- genometranslator::tidy_genome(data = genome)
```

The former `tidy_vcf()` wrapper has been removed. Its import
responsibilities belong to `read_genome()` and `read_vcf()`, while
filtering belongs in `radr`. The historical implementation remains
available in `radiator` for legacy workflows.

## Installation

`genometranslator` uses CRAN, Bioconductor, and GitHub dependencies.
Starting from a basic R installation, install the required dependencies
with:

``` r
install.packages(c("BiocManager", "remotes"))

BiocManager::install(c(
  "gdsfmt",
  "Rsamtools",
  "SeqArray"
))

remotes::install_github("thierrygosselin/tgbase")
remotes::install_github("thierrygosselin/genometranslator")
```

### Optional R packages

Install only the packages needed by your workflow:

``` r
# Interested in genind and/or genlight requires: 
install.packages("adegenet")
```

See the documentation for the relevant `read_*` or `write_*` function
for the exact dependency and data requirements of that format.

VCF imports retain caller and workflow provenance in GDS metadata and
preserve common structural-variant annotations when they are present.
See the [get-started vignette](articles/using_genometranslator.html) for
the current scope and the roadmap toward a dedicated structural-variant
representation.

Run the dependency diagnostic after installation:

``` r
genometranslator::genometranslator_dependencies()
```

### Optional bcftools executable

Some VCF preparation and indexing operations use `bcftools`. It is not
an R package and cannot be installed with `install.packages()`.

If Conda is not already available, install
[Miniforge](https://github.com/conda-forge/miniforge), restart the
terminal, and confirm that Conda is working:

``` bash
conda --version
```

Create the shared `genomics` environment when it does not already exist:

``` bash
conda create --name genomics --channel conda-forge --channel bioconda bcftools
```

For an existing `genomics` environment, install or update `bcftools`
with:

``` bash
conda activate genomics
conda install --channel conda-forge --channel bioconda bcftools
```

Confirm the installation before starting R or RStudio:

``` bash
conda activate genomics
bcftools --version
which bcftools
```

Start R or RStudio from the activated environment. Inside R, verify that
the executable is visible:

``` r
Sys.which("bcftools")
genometranslator::genometranslator_dependencies()
```

An empty result from `Sys.which("bcftools")` means that the current R
session cannot see the Conda environment. Close R/RStudio, activate
`genomics`, and start it again from that terminal session.

## Citation

To obtain the canonical citation for the installed package version, use:

``` r
citation("genometranslator")
```

When reporting an analysis, include at least the package version:

``` r
packageVersion("genometranslator")
```

For a development version, reproducibility is improved by also recording
the Git commit used. The access date is useful additional context,
especially while the package is under active development, but it should
not replace the version or commit identifier.

A citation can be described in this form until a dedicated publication
or DOI is available:

> Gosselin, T. (2026). *genometranslator: Read, standardize and
> translate genomic data*. R package version 0.0.0.9000.
> <https://github.com/thierrygosselin/genometranslator>. Accessed
> 2026-09-03.

## Website and support

Documentation and articles are available at
<https://thierrygosselin.github.io/genometranslator/>.

Report problems or request features through the [GitHub issue
tracker](https://github.com/thierrygosselin/genometranslator/issues).
