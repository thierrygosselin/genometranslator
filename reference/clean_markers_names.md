# Clean marker's names for genomic workflows

Function to clean marker's name of weird separators that interfere with
some packages or codes. `/`, `:`, `-` and `.` are changed to an
underscore `_`. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

## Usage

``` r
clean_markers_names(x)
```

## Arguments

- x:

  (character string) Markers character string.
