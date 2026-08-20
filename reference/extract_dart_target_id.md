# Extract [DArT](http://www.diversityarrays.com) TARGET_ID

Extracts DArT `TARGET_ID` from a DArT file to help build a STRATA file.
Optionally extracts the DArT metadata block.

DArT metadata are not fully consistent across files (fields present and
ordering may vary and are not always documented). When metadata are
requested, the function uses a conservative strategy:

- `TARGET_ID` is always identified (stable DArT position)

- `PLATE_WELL`, `WELL_ROW`, `WELL_COL` are detected when possible

- remaining metadata are named `DART_METADATA_*`

Possible DArT metadata fields may include:

- `DART_NUMBER`

- `DART_PLATE_BARCODE`

- `SAMPLE_COMMENTS`

- `WELL_ROW`

- `WELL_COL`

- `CLIENT_PLATE_BARCODE`

- `CLIENT_ID`

If known, metadata column names can be supplied using
`metadata.colnames`.

## Usage

``` r
extract_dart_target_id(
  data,
  write = TRUE,
  metadata = FALSE,
  metadata.colnames = NULL
)
```

## Arguments

- data:

  (file) 6 files formats used by DArT are recognized by radr. Don't
  modify the DArT file, to do this, use the `strata` file/argument
  below. The function can import files ending with `.csv` or `.tsv`.

  1.  **1row**: Genotypes are in 1 row and coded (0, 1, 2, -).
      `0 for 2 reference alleles REF/REF`,
      `1 for 2 alternate alleles ALT/ALT`, `2 for heterozygote REF/ALT`,
      `- for missing`.

  2.  **2rows**: No genotypes. It's absence/presence, 0/1, of the REF
      and ALT alleles. Sometimes called binary format.

  3.  **counts**: No genotypes, It's counts/read depth for the REF and
      ALT alleles. Sometimes just called count data. This should be the
      preferred file format, because DArT output the coverage (read
      depth for each genotypes).

  4.  **silico.dart**: SilicoDArT data. No genotypes, no REF or ALT
      alleles. It's a file coded as absence/presence, 0/1, for the
      presence of sequence in the clone id.

  5.  **silico.dart.counts**: SilicoDArT data. No genotypes, no REF or
      ALT alleles. It's a file coded as absence/presence, with counts
      for the presence of sequence in the clone id.

  6.  **dart.vcf**: For DArT VCFs, please use
      [`read_vcf`](https://thierrygosselin.github.io/genometranslator/reference/read_vcf.md).

  If you encounter a problem, sent me your data so that I can update the
  function.

- write:

  (logical) Default: `TRUE`. Write extracted TARGET_ID table.

  Default: `write = TRUE`.

- metadata:

  (logical) Default: `FALSE`. Extract metadata when present.

  Default: `metadata = FALSE`.

- metadata.colnames:

  (character, optional) Metadata column names excluding `TARGET_ID`.
  Length must match number of metadata columns. With the default,
  (automatic conservative naming).

  Default: `metadata.colnames = NULL`.

## Value

A tidy dataframe with `TARGET_ID` and optional metadata. TARGET_ID are
cleaned using
[`clean_ind_names()`](https://thierrygosselin.github.io/genometranslator/reference/clean_ind_names.md):
spaces and commas removed, `_` and `:` replaced by `-`, converted to
upper case.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>

## Examples

``` r
if (FALSE) { # \dontrun{

# Extract TARGET_ID only
genometranslator::extract_dart_target_id("mt.dart.file.csv")

# Extract metadata with automatic naming
genometranslator::extract_dart_target_id(
  data = "mt.dart.file.csv",
  metadata = TRUE
)

# Extract metadata with user naming
genometranslator::extract_dart_target_id(
  data = "mt.dart.file.csv",
  metadata = TRUE,
  metadata.colnames = c(
    "DART_NUMBER",
    "DART_PLATE_BARCODE",
    "CLIENT_PLATE_BARCODE",
    "WELL_ROW",
    "WELL_COL",
    "SAMPLE_COMMENTS"
  )
)

# Build STRATA file
strata <- genometranslator::extract_dart_target_id("mt.dart.file.csv") %>%
  dplyr::mutate(
    INDIVIDUALS = TARGET_ID,
    STRATA = "POP1"
  )

} # }
```
