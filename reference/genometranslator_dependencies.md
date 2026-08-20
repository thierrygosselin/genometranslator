# Check genometranslator dependencies

Reports required and optional R packages and checks whether the optional
`bcftools` executable is available. This function is intentionally
diagnostic: it does not modify the R library or a Conda environment.

## Usage

``` r
genometranslator_dependencies(verbose = TRUE)
```

## Arguments

- verbose:

  Logical. Print installation guidance for missing components. Default:
  `verbose = TRUE`.

## Value

A tibble with component, source, requirement level, and availability.
