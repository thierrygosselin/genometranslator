# Join stratification metadata to a VCF (population-aware VCF)

Include stratification metadata, e.g. population-level information, to
the `FORMAT` field of a VCF file.

## Usage

``` r
vcf_strata(data, strata, filename = NULL)
```

## Arguments

- data:

  A VCF file

- strata:

  (optional) A tab delimited file at least 2 columns with header:
  `INDIVIDUALS` and `STRATA`. The `STRATA` and any other columns can be
  any hierarchical grouping. To create a strata file see
  [`individuals2strata`](https://thierrygosselin.github.io/genometranslator/reference/individuals2strata.md).

- filename:

  (optional) The file name for the modifed VCF, written to the working
  directory. Default: `filename = NULL` will make a custom filename with
  data and time. Default: `filename = NULL`.

## Value

A VCF file in the working directory with new `FORMAT` field(s)
correponding to the strata column(s).

## References

Danecek P, Auton A, Abecasis G et al. (2011) The variant call format and
VCFtools. Bioinformatics, 27, 2156-2158.

## See also

[VCF web page](https://vcftools.github.io)

[https://vcftools.github.io/specs.html](https://thierrygosselin.github.io/genometranslator/reference/VCF%20specification%20page)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
