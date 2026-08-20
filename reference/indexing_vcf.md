# Check and ensure that a VCF is bgzipped and indexed

Internal helpers used by
[`read_vcf()`](https://thierrygosselin.github.io/genometranslator/reference/read_vcf.md)
to work with VCF files that are suitable for parallel reading with
SeqArray, i.e. bgzip- compressed and Tabix-indexed.

`detect_indexing()` checks whether a VCF file:

- is bgzip-compressed (`.vcf.gz`);

- has an existing Tabix index (`.tbi` or `.csi`);

- and that the index can be opened without error.

`indexing_vcf()` ensures that a VCF file is bgzip-compressed and indexed
with Tabix. If the file is not `.vcf.gz`, it is compressed with
[`Rsamtools::bgzip()`](https://rdrr.io/pkg/Rsamtools/man/zip.html); if
no index is found, it is created via
[`Rsamtools::indexTabix()`](https://rdrr.io/pkg/Rsamtools/man/indexTabix.html).

The typical workflow is to first call `detect_indexing()` and, if it
returns `TRUE`, to run `indexing_vcf()` to make the VCF compatible with
parallel chunked reading.

## Usage

``` r
detect_indexing(vcf)

indexing_vcf(vcf, bcftools.path = "bcftools", verbose = TRUE)
```

## Arguments

- vcf:

  (character) Path to a VCF file. Can be plain text (`.vcf`) or bgzipped
  (`.vcf.gz`).

- bcftools.path:

  Default: `bcftools.path = "bcftools"`.

- verbose:

  Default: `verbose = TRUE`.

## Value

- `detect_indexing()`: (logical) `TRUE` if the VCF needs compression
  and/or indexing (i.e. is not ready for parallel reading), `FALSE` if
  it is already bgzipped and indexed with a readable Tabix index.

- `indexing_vcf()`: (character) The path to the bgzipped and
  Tabix-indexed VCF.

## Details

These functions are designed for internal use and are not exported.
`indexing_vcf()` emits user-friendly progress messages via cli when
compression or indexing is required, but remains quiet when files are
already compliant.

## Examples

``` r
if (FALSE) { # \dontrun{
if (detect_indexing(vcf)) {
  vcf <- indexing_vcf(vcf)
}
} # }
```
