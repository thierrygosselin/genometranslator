# Read VCF files and write a radiator GDS file

Read a VCF file and convert it to a
[SeqArray](https://github.com/zhengxwen/SeqArray) GDS object
(`SeqVarGDSClass`; Zheng et al. 2017), with a `/metadata` node
containing genome-specific metadata.

The function has an "advanced" mode (via `...`) that allows several
VCF-specific clean-ups.

For users who want a fast and robust VCF-to-GDS import.

## Usage

``` r
read_vcf(
  data,
  strata = NULL,
  filename = NULL,
  vcf.stats = FALSE,
  parallel.core = parallel::detectCores() - 1,
  verbose = TRUE,
  ...
)
```

## Arguments

- data:

  (character) Path to a VCF file (optionally bgzipped, `*.vcf` or
  `*.vcf.gz`). Markers can be biallelic SNPs or haplotypes.

- strata:

  (optional) Strata definition, passed to
  [`read_strata`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md).
  Can be a path to a strata file or an object. Default: `strata = NULL`.

- filename:

  (optional, character) Base name of the GDS file to generate. Radiator
  will append `.gds.rad` to the filename. If the chosen filename already
  exists in `path.folder`, a timestamped default name is used instead.
  Default: `filename = NULL`.

- vcf.stats:

  (logical, optional) Generate basic statistics for individuals and
  markers (missingness, coverage, etc.) and write them to disk.
  Computational cost can be high for very large unfiltered VCF Default:
  `vcf.stats = FALSE`.

- parallel.core:

  (integer, optional) Approximate number of cores to use where
  parallelism is supported (e.g.
  [`SeqArray::seqApply`](https://rdrr.io/pkg/SeqArray/man/seqApply.html),
  some filter helpers). Default:
  `parallel.core = parallel::detectCores() - 1`.

- verbose:

  (logical, optional) When `TRUE`, the function prints progress messages
  and a summary. Default: `verbose = TRUE`.

## Value

A `SeqVarGDSClass` object (GDS) with a `/metadata` node containing
per-marker and per-individual metadata and a record of all filters
applied.

## Details

Typical performance (rough order of magnitude):

- a 35 GB VCF with ~4M SNPs: \\\sim\\7 minutes with 8 CPU;

- a 21 GB VCF with ~2M SNPs: \\\sim\\4 minutes with 8 CPU.

The resulting GDS file can be reopened almost instantly in a later R
session with
[`genometranslator::read_genome()`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.md).
So it's worth waiting.

## Dependencies

Required package dependencies are declared in DESCRIPTION and installed
with genometranslator. Run
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect core packages, optional packages, and external executables.

VCF import uses the declared Bioconductor packages SeqArray, gdsfmt, and
Rsamtools. VCF preparation and indexing may also require the optional
`bcftools` executable. See the installation section of the package
README and check visibility with `Sys.which("bcftools")`.

## Provenance and structural variants

The GDS records a VCF provenance table containing the original file
path, size, modification time, MD5 checksum, import time, inferred
caller, exact `##source` value, imported INFO and FORMAT fields, and
available reference, contig, command, and version header lines. This
information helps distinguish biological regional signals from caller,
mapper, reference-build, batch, or processing effects.

Common structural-variant INFO fields, including `SVTYPE`, `END`,
`SVLEN`, `CIPOS`, `CIEND`, `MATEID`, and `EVENT`, are retained when
present. Their presence does not make a biallelic dosage an adequate
representation of a structural variant. Structural variants and
candidate inversion-associated haploblocks should be annotated and
interpreted explicitly. A local genomic signal should be called a
candidate or putative inversion only until physical breakpoint, mapping,
assembly, linkage, or cytogenetic evidence supports it.

## VCF file format behaviour

**PLINK:**

- `LOCUS` is filled with an integer based on `CHROM`
  (`as.integer(factor(x = CHROM))`);

- `COL` is set to `1L` (no within-read position available).

**ipyrad:**

- the pattern `"locus_"` is stripped from `CHROM`;

- `COL` is set to `POS`.

**GATK / platypus / FreeBayes / samtools:**

- if the VCF `ID` column is `.`, it is replaced with the position
  (`POS`);

- short-read locus identity is assumed to be encoded in `CHROM` + `POS`.

**Stacks:**

- *de novo*: `CHROM` is typically "1"; `LOCUS` corresponds to "CHROM" in
  the Stacks VCF; `COL` is `POS - 1`;

- *reference*: `ID` is split into `LOCUS`, `COL`, `STRANDS`.

**DArT VCFs:**

- `CHROM == "."` is replaced by `"denovo"`;

- missing `POS` (`NA`) are set to `50`;

- `COL` is extracted from `LOCUS` (read position);

- `LOCUS` is the first group of digits, then joined with `POS` using
  `"_"`, and `POS` is replaced by `COL`.

## Advanced mode (`...`)

The `...` provides import and output controls that are intentionally
kept outside the basic interface:

- `blacklist.id`, `pop.select`, `pop.levels`, `pop.labels`;

- `filter.strands` - handle duplicate SNPs on opposite strands;

- `markers.info`, `vcf.metadata`;

- `path.folder`, `random.seed`, `subsample.markers.stats`;

- `filter.haplotype.format`, `parameters`, `internal`.

## References

Zheng X, Gogarten S, Lawrence M, Stilp A, Conomos M, Weir BS, Laurie C,
Levine D (2017). SeqArray – A storage-efficient high-performance data
format for WGS variant calls. *Bioinformatics*.

Danecek P, Auton A, Abecasis G et al. (2011) The variant call format and
VCFtools. *Bioinformatics* 27:2156-2158.

## See also

[`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md),
[`tidy_vcf`](https://thierrygosselin.github.io/genometranslator/reference/tidy_vcf.md)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{
# Simple import, no strata, defaults:
gds <- genometranslator::read_vcf(data = "populations.snps.vcf")

# With strata and a few filters:
gds <- genometranslator::read_vcf(
  data                      = "populations.snps.vcf",
  strata                    = "strata_salamander.tsv",
  path.folder               = "salamander",
  filter.strands             = "blacklist",
  verbose                    = TRUE
)

# Later, in a new R session, reopen the GDS:
gds <- genometranslator::read_genome(data = "radiator_20200911@0748.gds")
} # }
```
