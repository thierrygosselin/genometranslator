# detect_gt

Detect the genotype format used in the data set.

## Usage

``` r
detect_gt(
  x,
  gt.format = c("GT", "ALT_DOSAGE", "GT_VCF", "GT_VCF_NUC"),
  keep.one = TRUE,
  favorite = "ALT_DOSAGE"
)
```

## Arguments

- x:

  The data

- gt.format:

  (character) Default:
  `gt.format = c("GT", "ALT_DOSAGE", "GT_VCF", "GT_VCF_NUC")`.

- keep.one:

  (logical) Will return only one format if `keep.one = TRUE`. Default:
  `keep.one = TRUE`.

- favorite:

  If more than one format is present and `keep.one = TRUE`, the favorite
  will be returned, if present. Otherwise, the first format in
  `gt.format` is returned. Default: `favorite = "ALT_DOSAGE"`.
