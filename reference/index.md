# Package index

## Core workflow

Read, inspect, and translate genomic data.

- [`read_genome()`](https://thierrygosselin.github.io/genometranslator/reference/read_genome.md)
  : Read genomic data
- [`write_genome()`](https://thierrygosselin.github.io/genometranslator/reference/write_genome.md)
  : Write genomic data
- [`genome_translator()`](https://thierrygosselin.github.io/genometranslator/reference/genome_translator.md)
  : Translate genomic data between formats
- [`tidy_genome()`](https://thierrygosselin.github.io/genometranslator/reference/tidy_genome.md)
  : Convert a GDS genome to a tidy table
- [`genome_info()`](https://thierrygosselin.github.io/genometranslator/reference/genome_info.md)
  : Summarise genomic data dimensions
- [`genometranslator_dependencies()`](https://thierrygosselin.github.io/genometranslator/reference/genometranslator_dependencies.md)
  : Check genometranslator dependencies

## Sample metadata and strata

Prepare, validate, summarize, and modify sample groupings.

- [`read_strata()`](https://thierrygosselin.github.io/genometranslator/reference/read_strata.md)
  : read strata
- [`summary_strata()`](https://thierrygosselin.github.io/genometranslator/reference/summary_strata.md)
  : Summary of strata
- [`individuals2strata()`](https://thierrygosselin.github.io/genometranslator/reference/individuals2strata.md)
  : Create a strata file from a list of individuals
- [`generate_strata()`](https://thierrygosselin.github.io/genometranslator/reference/generate_strata.md)
  : Generate strata object from the data
- [`join_strata()`](https://thierrygosselin.github.io/genometranslator/reference/join_strata.md)
  : Join the strata with the data
- [`change_pop_names()`](https://thierrygosselin.github.io/genometranslator/reference/change_pop_names.md)
  : Transform into a factor the STRATA column, change names and reorder
  the levels
- [`clean_ind_names()`](https://thierrygosselin.github.io/genometranslator/reference/clean_ind_names.md)
  : Clean individual's names for genomic workflows
- [`clean_pop_names()`](https://thierrygosselin.github.io/genometranslator/reference/clean_pop_names.md)
  : Clean population's names for genomic workflows
- [`check_pop_levels()`](https://thierrygosselin.github.io/genometranslator/reference/check_pop_levels.md)
  : Check the use of pop.levels, pop.labels and pop.select arguments.
- [`vcf_strata()`](https://thierrygosselin.github.io/genometranslator/reference/vcf_strata.md)
  : Join stratification metadata to a VCF (population-aware VCF)

## Read genomic data

Format-specific readers offering additional control.

- [`read_dart()`](https://thierrygosselin.github.io/genometranslator/reference/read_dart.md)
  :

  Read and tidy [DArT](http://www.diversityarrays.com) output files.

- [`read_fstat()`](https://thierrygosselin.github.io/genometranslator/reference/read_fstat.md)
  : Read an FSTAT file into a tidy or wide data frame

- [`read_genepop()`](https://thierrygosselin.github.io/genometranslator/reference/read_genepop.md)
  : Read a Genepop file into a tidy or wide data frame

- [`read_genind()`](https://thierrygosselin.github.io/genometranslator/reference/read_genind.md)
  : Read a genind object to a GDS or tidy dataframe

- [`read_genlight()`](https://thierrygosselin.github.io/genometranslator/reference/read_genlight.md)
  : Read a genlight object into a tidy data frame and/or GDS object/file

- [`read_gtypes()`](https://thierrygosselin.github.io/genometranslator/reference/read_gtypes.md)
  : Read a gtypes object into a tidy data frame

- [`read_plink()`](https://thierrygosselin.github.io/genometranslator/reference/read_plink.md)
  : Reads PLINK tped and bed files

- [`read_vcf()`](https://thierrygosselin.github.io/genometranslator/reference/read_vcf.md)
  : Read VCF files and write a radiator GDS file

## Write genomic data

Export genomic data without silently filtering samples or markers.

- [`write_arlequin()`](https://thierrygosselin.github.io/genometranslator/reference/write_arlequin.md)
  : Write an arlequin file from a tidy data frame

- [`write_bayescan()`](https://thierrygosselin.github.io/genometranslator/reference/write_bayescan.md)
  :

  Write a [BayeScan](http://cmpg.unibe.ch/software/BayeScan/) file from
  a tidy data frame

- [`write_betadiv()`](https://thierrygosselin.github.io/genometranslator/reference/write_betadiv.md)
  : Write a betadiv file from a tidy data frame

- [`write_colony()`](https://thierrygosselin.github.io/genometranslator/reference/write_colony.md)
  :

  Write a `COLONY` input file

- [`write_dadi()`](https://thierrygosselin.github.io/genometranslator/reference/write_dadi.md)
  :

  Write a `dadi` SNP input file from a tidy data frame.

- [`write_faststructure()`](https://thierrygosselin.github.io/genometranslator/reference/write_faststructure.md)
  : Write a faststructure file from a tidy data frame

- [`write_fineradstructure()`](https://thierrygosselin.github.io/genometranslator/reference/write_fineradstructure.md)
  : Write a fineRADstructure file from a tidy data frame

- [`write_gds()`](https://thierrygosselin.github.io/genometranslator/reference/write_gds.md)
  : Write a GDS object from a tidy data frame

- [`write_genepop()`](https://thierrygosselin.github.io/genometranslator/reference/write_genepop.md)
  : Write a genepop file

- [`write_genepopedit()`](https://thierrygosselin.github.io/genometranslator/reference/write_genepopedit.md)
  : Write a genepopedit flatten object

- [`write_genind()`](https://thierrygosselin.github.io/genometranslator/reference/write_genind.md)
  : Write a genind object from a tidy data frame or GDS file or object.

- [`write_genlight()`](https://thierrygosselin.github.io/genometranslator/reference/write_genlight.md)
  :

  Write a `genlight` object from: a tidy data frame, GDS file or object.

- [`write_genome()`](https://thierrygosselin.github.io/genometranslator/reference/write_genome.md)
  : Write genomic data

- [`write_gsi_sim()`](https://thierrygosselin.github.io/genometranslator/reference/write_gsi_sim.md)
  : Write a gsi_sim file from a data frame (wide or long/tidy).

- [`write_gtypes()`](https://thierrygosselin.github.io/genometranslator/reference/write_gtypes.md)
  : Write a strataG gtypes object from GDS or tidy data

- [`write_hapmap()`](https://thierrygosselin.github.io/genometranslator/reference/write_hapmap.md)
  : Write a HapMap file from a tidy data frame

- [`write_hierfstat()`](https://thierrygosselin.github.io/genometranslator/reference/write_hierfstat.md)
  : Write a hierfstat file from a tidy data frame

- [`write_hzar()`](https://thierrygosselin.github.io/genometranslator/reference/write_hzar.md)
  : Write a HZAR file from a tidy data frame.

- [`write_ldna()`](https://thierrygosselin.github.io/genometranslator/reference/write_ldna.md)
  : Write a LDna object from a tidy data frame

- [`write_maverick()`](https://thierrygosselin.github.io/genometranslator/reference/write_maverick.md)
  : Write a maverick file from a tidy data frame

- [`write_pcadapt()`](https://thierrygosselin.github.io/genometranslator/reference/write_pcadapt.md)
  :

  Write a [pcadapt](https://github.com/bcm-uga/pcadapt) file from a tidy
  data frame

- [`write_plink()`](https://thierrygosselin.github.io/genometranslator/reference/write_plink.md)
  : Write a plink tped/tfam file from a tidy data frame

- [`write_related()`](https://thierrygosselin.github.io/genometranslator/reference/write_related.md)
  : Write a related file from a tidy data frame

- [`write_rubias()`](https://thierrygosselin.github.io/genometranslator/reference/write_rubias.md)
  : Write a rubias object

- [`write_snprelate()`](https://thierrygosselin.github.io/genometranslator/reference/write_snprelate.md)
  : Write a SNPRelate object from a tidy data frame

- [`write_stockr()`](https://thierrygosselin.github.io/genometranslator/reference/write_stockr.md)
  : Write a stockR dataset from a tidy data frame or GDS file or object.

- [`write_structure()`](https://thierrygosselin.github.io/genometranslator/reference/write_structure.md)
  : Write a structure file from a tidy data frame

- [`write_vcf()`](https://thierrygosselin.github.io/genometranslator/reference/write_vcf.md)
  : Write a vcf file from a tidy data frame

## Detect formats and encoding

- [`detect_genomic_format()`](https://thierrygosselin.github.io/genometranslator/reference/detect_genomic_format.md)
  : Used internally in radiator to detect the file format
- [`detect_dart_format()`](https://thierrygosselin.github.io/genometranslator/reference/detect_dart_format.md)
  : detect_dart_format
- [`detect_gt()`](https://thierrygosselin.github.io/genometranslator/reference/detect_gt.md)
  : detect_gt
- [`detect_biallelic_markers()`](https://thierrygosselin.github.io/genometranslator/reference/detect_biallelic_markers.md)
  : Detect biallelic data
- [`detect_indexing()`](https://thierrygosselin.github.io/genometranslator/reference/indexing_vcf.md)
  [`indexing_vcf()`](https://thierrygosselin.github.io/genometranslator/reference/indexing_vcf.md)
  : Check and ensure that a VCF is bgzipped and indexed

## Work with GDS

- [`genome_gds()`](https://thierrygosselin.github.io/genometranslator/reference/genome_gds.md)
  : Genome GDS constructor
- [`genome_gds_skeleton()`](https://thierrygosselin.github.io/genometranslator/reference/genome_gds_skeleton.md)
  : genome_gds_skeleton
- [`summary_gds()`](https://thierrygosselin.github.io/genometranslator/reference/summary_gds.md)
  : summary_gds
- [`update_genome_gds()`](https://thierrygosselin.github.io/genometranslator/reference/update_genome_gds.md)
  : update_genome_gds
- [`upgrade_genome_gds()`](https://thierrygosselin.github.io/genometranslator/reference/upgrade_genome_gds.md)
  : Upgrade a legacy radiator GDS file
- [`sync_gds()`](https://thierrygosselin.github.io/genometranslator/reference/sync_gds.md)
  : sync_gds
- [`tidy2gds()`](https://thierrygosselin.github.io/genometranslator/reference/tidy2gds.md)
  : tidy2gds
- [`parse_gds_metadata()`](https://thierrygosselin.github.io/genometranslator/reference/parse_gds_metadata.md)
  : parse_gds_metadata
- [`genome_parameters()`](https://thierrygosselin.github.io/genometranslator/reference/genome_parameters.md)
  : Track changes to genomic data
- [`list_filters()`](https://thierrygosselin.github.io/genometranslator/reference/list_filters.md)
  : List current active filters (individuals and markers) in radiator
  GDS object.
- [`reset_filters()`](https://thierrygosselin.github.io/genometranslator/reference/reset_filters.md)
  : Reset filters (individuals and markers) in radiator GDS object.

## DArT utilities

- [`extract_dart_target_id()`](https://thierrygosselin.github.io/genometranslator/reference/extract_dart_target_id.md)
  :

  Extract [DArT](http://www.diversityarrays.com) TARGET_ID

- [`extract_dart_markers_metadata()`](https://thierrygosselin.github.io/genometranslator/reference/extract_dart_markers_metadata.md)
  : extract_dart_markers_metadata

- [`tidy_dart_metadata()`](https://thierrygosselin.github.io/genometranslator/reference/tidy_dart_metadata.md)
  :

  Import and tidy [DArT](http://www.diversityarrays.com) metadata.

## VCF and bcftools utilities

- [`bcftools_require()`](https://thierrygosselin.github.io/genometranslator/reference/bcftools_require.md)
  : Check that bcftools is available
- [`bcftools_exec()`](https://thierrygosselin.github.io/genometranslator/reference/bcftools_exec.md)
  : Run a bcftools command and log stderr
- [`detect_indexing()`](https://thierrygosselin.github.io/genometranslator/reference/indexing_vcf.md)
  [`indexing_vcf()`](https://thierrygosselin.github.io/genometranslator/reference/indexing_vcf.md)
  : Check and ensure that a VCF is bgzipped and indexed
- [`check_header_source_vcf()`](https://thierrygosselin.github.io/genometranslator/reference/check_header_source_vcf.md)
  : Check a VCF header and detect its source (caller)
- [`extract_individuals_vcf()`](https://thierrygosselin.github.io/genometranslator/reference/extract_individuals_vcf.md)
  : Extract individuals from vcf file

## Genotype and marker utilities

- [`extract_genotypes_metadata()`](https://thierrygosselin.github.io/genometranslator/reference/extract_genotypes_metadata.md)
  : extract_genotypes_metadata
- [`extract_individuals_metadata()`](https://thierrygosselin.github.io/genometranslator/reference/extract_individuals_metadata.md)
  : extract_individuals_metadata
- [`extract_markers_metadata()`](https://thierrygosselin.github.io/genometranslator/reference/extract_markers_metadata.md)
  : extract_markers_metadata

## Whitelists and blacklists

Import lists used to retain or exclude samples and markers.

- [`read_blacklist_id()`](https://thierrygosselin.github.io/genometranslator/reference/read_blacklist_id.md)
  : read_blacklist_id
- [`read_whitelist()`](https://thierrygosselin.github.io/genometranslator/reference/read_whitelist.md)
  : Read a marker whitelist
