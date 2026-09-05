
# genometranslator <a href="https://thierrygosselin.github.io/genometranslator/"><img src="man/figures/logo.png" align="right" height="160" /></a>

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![packageversion](https://img.shields.io/badge/Package%20version-0.0.0.9000-orange.svg)](commits/master)
[![Last-changedate](https://img.shields.io/badge/last%20change-2026--09--05-brightgreen.svg)](/commits/master)
<!-- badges: end -->

`genometranslator` reads, standardizes, and writes individual genomic
data. The easiest interface is `read_genome()` and `write_genome()`. Use
the format-specific `read_*()` or `write_*()` functions when more
control is needed.

The genomic datasets remain in a file-backed GDS rather than being
expanded into an R data frame unless that conversion is requested
explicitly.

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
  vcf.stats = TRUE,
  parallel.core = 20L,
  all.sites = TRUE
)
```

Use `tidy_genome()` only when an in-memory table is genuinely needed:

``` r
genotypes <- genometranslator::tidy_genome(data = genome)
```

## Prepare an OutFLANK scan

`write_outflank()` calculates corrected and uncorrected population FST
components directly from diploid biallelic GDS calls, reading loci in
blocks. It writes the OutFLANK statistics table and sample, population
and locus audits. Filenames include a `YYYYMMDD@HHMM` timestamp;
existing files are never replaced. Population metadata are matched by
sample ID. The default stops on insufficient per-population calls or
monomorphic loci; explicit exclusions are audited.

``` r
exported <- genometranslator::write_outflank(
  data = "study.gds", strata = "samples.tsv", filename = "study"
)
```

See `?write_outflank` for assumptions, references and the downstream
example. The writer does not run a selection scan or fix downstream
OutFLANK issues.

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

Run the dependency diagnostic after installation:

``` r
genometranslator::genometranslator_dependencies()
```

See the documentation for the relevant `read_*` or `write_*` function
for the exact dependency and data requirements of that format.

More information in the [get-started
vignette](articles/using_genometranslator.html).

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
> 2026-09-05.

## Website and support

### Choosing a converter

Read [Choosing a genomic data
converter](https://thierrygosselin.github.io/genometranslator/articles/converter_comparison.html)
for tested conversion paths, known limitations, and genometranslator’s
own safeguards. Format support alone does not guarantee faithful
conversion.

### Requesting additional genomic formats

Need a genomic format that genometranslator does not support? Submit a
[format request on
GitHub](https://github.com/thierrygosselin/genometranslator/issues).
Include the format specification, the intended workflow, and a small
synthetic or non-sensitive example. Requests are evaluated for
scientific usefulness, reliable translation, and maintenance
requirements; support is not guaranteed.

Dedicated **fastSTRUCTURE support will not be implemented**. See
[Unsupported formats and software-specific
exports](https://thierrygosselin.github.io/genometranslator/articles/unsupported_formats.html)
for the code-audit findings behind this decision. This vignette records
exclusions and their rationale as the package evolves. The
fineRADstructure writer is also retired; a request may prompt
reassessment of an updated implementation against the documented
findings. The related writer is also retired; the same vignette
documents the implementation concerns behind this support decision.

Documentation and articles are available at
<https://thierrygosselin.github.io/genometranslator/>.

Report problems or request features through the [GitHub issue
tracker](https://github.com/thierrygosselin/genometranslator/issues).
