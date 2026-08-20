# Tidy vcf file

The function allows to tidy a VCF file.

Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

It is highly recommended to use
[`radr::explore_genomes`](https://thierrygosselin.github.io/radr/reference/explore_genomes.html)
to reduce the number of markers. Advance options below are also
available to to manipulate and prune the dataset with blacklists and
whitelists along several other filtering options.

## Usage

``` r
tidy_vcf(
  data,
  strata = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1,
  verbose = FALSE,
  ...
)
```

## Arguments

- data:

  (VCF file, character string) The VCF SNPs are biallelic or haplotypes.
  To make the VCF population-ready, the argument `strata` is required.

- filename:

  (optional) The function uses
  [`read_parquet`](https://arrow.apache.org/docs/r/reference/read_parquet.html),
  to write the tidy data frame in the working directory. The file
  extension appended to the `filename` provided is `.arrow.parquet`.
  With default: `filename = NULL`, the tidy data frame is in the global
  environment only (i.e. not written in the working directory...).
  Default: `filename = NULL`.

- parallel.core:

  Number of workers available for parallel operations. Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  Logical indicating whether progress messages are emitted. Default:
  `verbose = FALSE`.

- ...:

  (optional) To pass further argument for fine-tuning the tidying (read
  below).

## Value

The output in your global environment is a tidy data frame, the GDS file
generated is in the working directory under the name given during
function execution.

## VCF file format

**PLINK:** radiator fills the `LOCUS` column of PLINK VCFs with a unique
integer based on the `CHROM` column (`as.integer(factor(x = CHROM))`).
The `COL` column is filled with 1L for lack of bettern info on this. Not
what you need ? Open an issue on GitHub for a request.

**ipyrad:** the pattern `locus_` in the `CHROM` column is removed and
used. The `COL` column is filled with the same value as `POS`.

**GATK:** Some VCF have an `ID` column filled with `.`, the LOCUS
information is all contained along the linkage group in the `CHROM`
column. To make it work with
[genometranslator](https://github.com/thierrygosselin/genometranslator),
the `ID` column is filled with the `POS` column info. GATK with a mix of
multi- and bi-allelic dataset won't generate VCF stats.

**platypus:** Some VCF files don't have an ID filed with values, here
the same thing is done as GATK VCF files above.

**freebayes:** Some VCF files don't have an ID filed with values, here
the same thing is done as GATK VCF files above.

**stacks:** with *de novo* approaches, the CHROM column is filled with
"1", the LOCUS column correspond to the CHROM section in stacks VCF and
the COL column is POS -1. With a reference genome, the ID column in
stacks VCF is separated into "LOCUS", "COL", "STRANDS".

*stacks problem:* current version as some intrinsic problem with missing
allele depth info, during the tidying process a message will highlight
the number of genotypes impacted by the problem. When possible, the
problem is corrected by adding the read depth info into the allele depth
field.

## Advance mode, using *dots-dots-dots ...*

The arguments below are not available using code completion (e.g. with
TAB), consequently any misspelling will generate an error or be ignored.

*dots-dots-dots ...* arguments names and values are reported and written
in the working directory.

**General arguments:**

1.  `path.folder`: to write ouput in a specific path (used internally in
    radiator). Default: `path.folder = getwd()`. If the supplied
    directory doesn't exist, it's created.

2.  `random.seed`: (integer, optional) For reproducibility, set an
    integer that will be used inside codes that uses randomness. With
    default, a random number is generated, printed and written in the
    appropriate directory. Random seed is recycled inside the function
    that will import the VCF file before tidying. Default:
    `random.seed = NULL`.

**tidying arguments/behavior:**

1.  `tidy.vcf:` (optional, logical) Default: `tidy.vcf = TRUE`. But you
    can always stop the process after the creation of the GDS file
    (equivalent of running
    [`read_vcf`](https://thierrygosselin.github.io/genometranslator/reference/read_vcf.md)).

2.  `tidy.check:` (optional, logical) Default: `tidy.check = TRUE`. By
    default, the number of markers just before tidying is checked.
    Tidying a VCF file with more than 20000 markers is sub-optimal:

    - a computer with lots of RAM is required

    - it's very slow to generate

    - it's very slow to run codes after

    - for most non model species this number of markers is not
      realistic...

    Consequently, the function execution is suspended and user are asked
    if they still want to continue with the tidying or stop and keep the
    GDS file/object.

    This behavior can be annoying, *if the user knows what he's doing*,
    to turn off use: `tidy.check = FALSE`.

3.  `calibrate.alleles: ` (optional, logical) Default:
    `calibrate.alleles = FALSE`. Documented in
    [`calibrate_alleles`](https://thierrygosselin.github.io/genometranslator/reference/calibrate_alleles.md).

4.  `vcf.stats: ` (optional, logical) Generates individuals and markers
    statistics helpful for filtering. These are very fast to generate
    and because computational cost is minimal, even for huge VCFs, the
    default is `vcf.stats = TRUE`.

5.  `vcf.metadata: ` (optional, logical or character string) With
    `vcf.metadata = FALSE`, only the genotypes are kept (GT field) in
    the tidy dataset. With `vcf.metadata = TRUE`, all the metadata
    contained in the `FORMAT` field will be kept in the tidy data file.
    radiator is currently keeping and cleaning these metadata:
    `"DP", "AD", "GL", "PL", "GQ", "HQ", "GOF", "NR", "NV", "CATG"`.
    e.g. you only want AD and PL, `vcf.metadata = c("AD", "PL")`. Need
    another metadata ? Submit a request on github... Default:
    `vcf.metadata = TRUE`.

**Filtering arguments:**

1.  `blacklist.id: ` (optional, character) Default
    (`blacklist.id = NULL`). Documented in
    [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

2.  `filter.strands`: (optional, character) Default
    (`filter.strands = "blacklist"`). documented in
    [`read_vcf`](https://thierrygosselin.github.io/genometranslator/reference/read_vcf.md).

3.  `whitelist.markers: `(optional, path) Default:
    `whitelist.markers = NULL`. Documented in
    [`radr::filter_whitelist`](https://thierrygosselin.github.io/radr/reference/filter_whitelist.html).

4.  `filter.individuals.missing`: (double) Default:
    `filter.individuals.missing = NULL`. Documented in
    [`radr::filter_individuals`](https://thierrygosselin.github.io/radr/reference/filter_individuals.html).

5.  `filter.monomorphic`: (logical) Default:
    `filter.monomorphic = TRUE`. Documented in
    [`radr::filter_monomorphic`](https://thierrygosselin.github.io/radr/reference/filter_monomorphic.html).
    Required package: `UpSetR`.

6.  `filter.common.markers`: (logical) Default:
    `filter.common.markers = TRUE`. Documented in
    [`radr::filter_common_markers`](https://thierrygosselin.github.io/radr/reference/filter_common_markers.html).
    Required package: `UpSetR`.

7.  `filter.ma`: (integer) Default: `filter.ma = NULL`. Documented in
    [`radr::filter_ma`](https://thierrygosselin.github.io/radr/reference/filter_ma.html).

8.  `filter.coverage`: (logical) Default: `filter.coverage = NULL`.
    Documented in
    [`radr::filter_coverage`](https://thierrygosselin.github.io/radr/reference/filter_coverage.html).

9.  `filter.genotyping`: (integer) Default: `filter.genotyping = NULL`.
    Documented in
    [`radr::filter_genotyping`](https://thierrygosselin.github.io/radr/reference/filter_genotyping.html).

10. `filter.snp.position.read: ` (optional, character, integer) Default:
    `filter.snp.position.read = NULL`. Documented in
    [`radr::filter_snp_position_read`](https://thierrygosselin.github.io/radr/reference/filter_snp_position_read.html).

11. `filter.snp.number: ` (optional, character, integer) Default:
    `filter.snp.number = NULL`. Documented in
    [`radr::filter_snp_number`](https://thierrygosselin.github.io/radr/reference/filter_snp_number.html).

12. `filter.short.ld`: (optional, character) Default:
    `filter.short.ld = NULL`. Documented in
    [`radr::filter_ld`](https://thierrygosselin.github.io/radr/reference/filter_ld.html).

13. `filter.long.ld`: (optional, character) Default:
    `filter.long.ld = NULL`. Documented in
    [`radr::filter_ld`](https://thierrygosselin.github.io/radr/reference/filter_ld.html).
    Required package: `SNPRelate`.

14. `long.ld.missing`: Documented in
    [`radr::filter_ld`](https://thierrygosselin.github.io/radr/reference/filter_ld.html).
    Default: `long.ld.missing = FALSE`.

15. `ld.method`: Documented in
    [`radr::filter_ld`](https://thierrygosselin.github.io/radr/reference/filter_ld.html).
    Default: `ld.method = "r2"`.

## References

Danecek P, Auton A, Abecasis G et al. (2011) The variant call format and
VCFtools. Bioinformatics, 27, 2156-2158.

## See also

[`read_vcf`](https://thierrygosselin.github.io/genometranslator/reference/read_vcf.md),
[`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md),
[`radr::explore_genomes`](https://thierrygosselin.github.io/radr/reference/explore_genomes.html)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# very basic with built-in defaults (not recommended):
prep.data <- genometranslator::tidy_vcf(data = "populations.snps.vcf")

# Using more arguments and filters (recommended):
tidy.data <- genometranslator::tidy_vcf(
    data = "populations.snps.vcf",
    strata = "strata_salamander.tsv",
    filter.individuals.missing = "outlier",
    filter.ma = 4,
    filter.genotyping = 0.1,
    filter.snp.position.read = "outliers",
    filter.short.ld = "mac",
    path.folder = "salamander/prep_data",
    verbose = TRUE)
} # }
```
