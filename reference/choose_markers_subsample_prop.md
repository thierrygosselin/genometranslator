# choose markers subsample prop

Choose the best subsampling prop based on different metrics

## Usage

``` r
choose_markers_subsample_prop(
  n.markers,
  n.ind,
  vector.size.limit = 2^31,
  max.mem.gb = 10,
  bytes.per.cell = 4,
  n.matrices = 1,
  small.vcf.cutoff = 2e+05,
  cap.prop = 0.3,
  allowed = c(1, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01),
  safety.margin = 0.9
)
```

## Arguments

- vector.size.limit:

  Default: `vector.size.limit = 2^31`.

- max.mem.gb:

  Default: `max.mem.gb = 10`.

- bytes.per.cell:

  Default: `bytes.per.cell = 4`.

- n.matrices:

  Default: `n.matrices = 1`.

- small.vcf.cutoff:

  Default: `small.vcf.cutoff = 2e+05`.

- cap.prop:

  Default: `cap.prop = 0.3`.

- allowed:

  Default: `allowed = c(1, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01)`.

- safety.margin:

  Default: `safety.margin = 0.9`.

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
