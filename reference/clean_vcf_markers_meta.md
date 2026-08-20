# title Clean and normalise markers metadata from a VCF import

title Clean and normalise markers metadata from a VCF import

## Usage

``` r
clean_vcf_markers_meta(
  markers.meta,
  gds,
  data.source,
  ref.genome,
  stacks.checks = FALSE,
  verbose = TRUE
)
```

## Arguments

- markers.meta:

  A tibble/data.frame of markers metadata containing at least
  `VARIANT_ID`, `CHROM`, `LOCUS`, `POS`.

- gds:

  A SeqArray GDS object with the VCF data already imported.

- data.source:

  (character) A short tag describing the originating pipeline (e.g.
  `"Stacks"`, `"PLINK"`, `"ipyrad"`, `"freeBayes"`). May be updated
  inside the function (e.g. to `"dart.vcf"` if the structure matches
  DArT VCF).

- ref.genome:

  (logical) Output of `detect_ref_genome()`, indicating whether the data
  were produced with a reference genome.

- stacks.checks:

  (logical) Whether to apply additional Stacks-specific adjustments when
  `ref.genome = FALSE`. Default: `stacks.checks = FALSE`.

- verbose:

  (logical) Display messages during cleaning. Default: `verbose = TRUE`.

## Value

A list with:

- `markers.meta`: cleaned tibble with harmonised fields;

- `data.source`: (possibly updated) data source tag;

- `detect.strand`: logical flag indicating whether strand info was
  detected from the LOCUS field.
