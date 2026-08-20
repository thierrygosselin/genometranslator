# Check a VCF header and detect its source (caller)

Inspect the VCF header (via
[`SeqArray::seqVCF_Header()`](https://rdrr.io/pkg/SeqArray/man/seqVCF_Header.html))
to:

- fix known FORMAT/INFO inconsistencies (Stacks,
  samtools/freebayes-style),

- detect the variant caller / pipeline (Stacks, freebayes, bcftools,
  GATK, etc.),

- flag whether the VCF comes from Stacks,

- suggest reasonable INFO/FORMAT subsets for import,

- return a cleaned header object for safe import into GDS.

The function preserves the original `##source=` line (if present) and
returns it as `vcf.source.raw`. Classification into a standardised
descriptor is returned as `data.source`.

## Usage

``` r
check_header_source_vcf(vcf, markers.info = NULL, vcf.metadata = NULL)
```

## Arguments

- vcf:

  (character) Path to a VCF file (optionally bgzipped).

- markers.info:

  (optional, character) Names of INFO fields to import. If `NULL`,
  caller-specific defaults are used when available; if those are also
  `NULL`, all INFO fields are imported. Any values not present in the
  header are silently dropped.

  Default: `markers.info = NULL`.

- vcf.metadata:

  (optional) Controls which FORMAT fields are imported:

  - `NULL`: use caller-specific defaults (if any), otherwise import all
    FORMAT fields;

  - logical: `TRUE` = import all FORMAT fields, `FALSE` = import only
    `GT`;

  - character: explicit list of FORMAT field IDs to import; `GT` is
    always added if missing. Any fields not present in the header are
    silently dropped.

  Default: `vcf.metadata = NULL`.

## Value

A list with:

- `data.source` – inferred caller (e.g. "bcftools", "freebayes", …);

- `vcf.source.raw` – exact raw `##source=` string or `NA`;

- `stacks.check` – whether this is a Stacks VCF (TRUE/FALSE);

- `check.header` – cleaned
  [`SeqArray::seqVCF_Header()`](https://rdrr.io/pkg/SeqArray/man/seqVCF_Header.html)
  object;

- `markers.info` – final INFO fields to import (or `NULL` = all);

- `overwrite.metadata` – final FORMAT fields to import (or `NULL` =
  all).
