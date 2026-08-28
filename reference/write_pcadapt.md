# Write a [pcadapt](https://github.com/bcm-uga/pcadapt) file from a tidy data frame

Write a [pcadapt](https://github.com/bcm-uga/pcadapt) file from a tidy
data frame. The data is biallelic. Used internally in
[genometranslator](https://github.com/thierrygosselin/genometranslator)
and might be of interest for users.

Prepare an appropriately filtered biallelic dataset before export.
Consider missingness, minor-allele frequency, linkage disequilibrium,
and the sampling design intended for the pcadapt analysis.

pcadapt candidates should be interpreted as loci unusually associated
with inferred population structure, not as confirmed targets of
selection. Use false-discovery-rate control on the p-values produced by
pcadapt, for example with the Bioconductor qvalue package, and state the
chosen q-value threshold. A threshold of 0.05 is common for candidate
discovery; 0.01 is a more conservative option when false positives are
costly.

## Usage

``` r
write_pcadapt(
  data,
  pop.select = NULL,
  filename = NULL,
  parallel.core = parallel::detectCores() - 1
)
```

## Arguments

- data:

  A tidy data frame object in the global environment or a tidy data
  frame in wide or long format in the working directory. *How to get a
  tidy data frame ?* Look into genometranslator
  [`tidy_genome`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md).

- pop.select:

  (optional, string) Selected list of populations for the analysis. e.g.
  `pop.select = c("QUE", "ONT")` to select `QUE` and `ONT` population
  samples (out of 20 pops). If `pop.labels` argument was used to rename
  the strata column, use the new names with `pop.select`. Default:
  `pop.select = NULL`.

- filename:

  (optional) The file name prefix for the pcadapt file written to the
  working directory. With default: `filename = NULL`, the date and time
  is appended to `radiator_pcadapt_`. Default: `filename = NULL`.

- parallel.core:

  Default: `parallel.core = parallel::detectCores() - 1`.

## Value

A pcadapt file is written in the working directory a genotype matrix
object is also generated in the global environment.

## Details

**Use a filtered dataset**:

1.  **Control linkage disequilibrium**: Reducing linkage before running
    genome scan is essential. At least start by removing SNPs on the
    same RADseq locus (short linkage disequilibrium).

2.  **Control minor alleles**: Too much Minor Alleles is just noise in
    the data. Filter using Count, Frequency or Depth.

3.  **Only use markers that are in common between strata**

4.  **Only use polymorphic markers**

**Use complementary genome scans**: No single genome-scan method is
robust to every demographic history, sampling design, or genetic
architecture. Compare pcadapt with methods based on different
assumptions when possible. For example,
[radr::run_bayescan()](https://thierrygosselin.github.io/radr/reference/run_bayescan.html)
implements a population-based FST outlier scan, whereas pcadapt models
individual genotypes through principal components without requiring
predefined populations. Environmental-association and haplotype-aware
approaches can add complementary evidence when the study design supports
them.

Published method comparisons must be interpreted within their evaluation
scenarios. The BayeScan false-discovery and admixture results reported
by Luu et al. (2017), for example, were obtained from simulations and
are not a universal ranking of the methods. Whole-genome analyses by
Meisner et al. (2021) also showed that pcadapt statistics can be
inflated under discrete population structure and when principal
components reflect sequencing or genotype-calling artefacts. Inspect the
PCA, QQ plots, genomic inflation, sample provenance, and technical
covariates before interpreting candidates.

Agreement among methods can increase confidence, but a strict
intersection is not automatically the best candidate set because the
methods detect different signals and have different failure modes.
Report method-specific results, examine sensitivity to filtering,
linkage, the number of principal components, missingness, and sampling,
and validate important candidates using independent data or
study-specific simulations.

**Candidate inversions and structural regions:** A large
inversion-associated or low-recombination haploblock can dominate a
principal component and consequently dominate a pcadapt scan. Compare
the complete genome, a collinear sensitivity dataset excluding candidate
regions, and an inversion-specific analysis. A local PCA signal may also
reflect a centromere, assembly or mapping problem, introgression,
population-specific missingness, or another structural variant. Because
genometranslator does not depend on radr, candidate screening and
genomic context are documented at [radr's inversion
vignette](https://thierrygosselin.github.io/radr/articles/detecting_inversions.html).

## Data filtering

This writer does not silently filter markers or individuals. It may
validate requirements imposed by the destination format and stop with an
informative error when the input is unsuitable. It is the user's
responsibility to filter and quality-control the data appropriately for
the intended analysis before generating the output. Use
[radr](https://thierrygosselin.github.io/radr/) or another suitable
workflow when filtering is required.

## Dependencies

Required package dependencies are declared in `DESCRIPTION` and are
installed with genometranslator. Any additional dependency needed only
for this format or option is identified in this help page. Use
[`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
to inspect the availability of core packages, optional packages, and
external executables.

## References

Luu, K., Bazin, E., & Blum, M. G. (2017). pcadapt: an R package to
perform genome scans for selection based on principal component
analysis. Molecular Ecology Resources, 17(1), 67-77.
[doi:10.1111/1755-0998.12592](https://doi.org/10.1111/1755-0998.12592)

Duforet-Frebourg, N., Luu, K., Laval, G., Bazin, E., & Blum, M. G.
(2015). Detecting genomic signatures of natural selection with principal
component analysis: application to the 1000 Genomes data. Molecular
biology and evolution, msv334.
[doi:10.1093/molbev/msv334](https://doi.org/10.1093/molbev/msv334)

de Villemereuil, P., Frichot, E., Bazin, E., Francois, O., & Gaggiotti,
O. E. (2014). Genome scan methods against more complex models: when and
how much should we trust them? Molecular Ecology, 23, 2006-2019.
[doi:10.1111/mec.12705](https://doi.org/10.1111/mec.12705)

Lotterhos, K. E., & Whitlock, M. C. (2014). Evaluation of demographic
history and neutral parameterization on the performance of FST outlier
tests. Molecular Ecology, 23, 2178-2192.
[doi:10.1111/mec.12725](https://doi.org/10.1111/mec.12725)

Meisner, J., Albrechtsen, A., & Hanghoj, K. (2021). Detecting selection
in low-coverage high-throughput sequencing data using principal
component analysis. BMC Bioinformatics, 22, 470.
[doi:10.1186/s12859-021-04375-2](https://doi.org/10.1186/s12859-021-04375-2)

## Author

Thierry Gosselin <thierrygosselin@icloud.com>
