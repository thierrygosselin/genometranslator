# Check that bcftools is available

Internal helper that checks whether `bcftools` can be executed. It runs
`bcftools.path --version` and errors if the command fails.

## Usage

``` r
bcftools_require(bcftools.path = "bcftools")
```

## Arguments

- bcftools.path:

  (character) Path or name of the `bcftools` executable. Default:
  `bcftools.path = "bcftools"`.
