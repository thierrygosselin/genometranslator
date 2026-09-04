#' Write Lep-MAP3
#'
#' @description Export a validated VCF and transposed pedigree for Lep-MAP3
#' ParentCall2, using genotype likelihoods retained in a SeqArray GDS.
#' @param data Open SeqArray GDS or existing GDS path.
#' @param pedigree Data frame or tab-separated file with FAMILY, INDIVIDUALS,
#' FATHER, MOTHER and SEX. IDs must match the active GDS samples exactly.
#' Parent IDs use 0 for unknown; SEX uses 1 (male), 2 (female), or 0 (unknown).
#' Optional PHENOTYPE defaults to 0 and is not used for linkage mapping.
#' @param likelihood Likelihood field: auto prefers PL, otherwise GL. An explicit
#' choice exports only that field alongside GT; no likelihoods are manufactured.
#' @param filename Output basename.
#' @param path.folder Output directory.
#' @param overwrite Replace existing export files.
#' @param verbose Display messages wrapped to 80 characters.
#' @return Invisibly, named paths to the VCF, pedigree, sample and marker mapping
#' tables, and a text file containing an example ParentCall2 command.
#' @details
#' \itemize{
#' \item Supports diploid, biallelic A/C/G/T SNPs. Unsupported variants are
#' rejected, not silently discarded. Active GDS filters and metadata whitelists
#' are respected; the caller's GDS selection is restored on success or error.
#' \item PL values must be non-negative integers; GL values must be finite
#' log10 likelihoods. Each observed vector contains three values in VCF order
#' (REF/REF, REF/ALT, ALT/ALT). Fully missing vectors remain missing; partial
#' vectors and markers without any observed likelihoods are rejected.
#' \item The minimal VCF retains GT and the selected likelihood field, not
#' unrelated annotations. PL/GL values are not converted to posterior calls.
#' \item Pedigrees are matched by sample ID, never by row order or STRATA.
#' This initial exporter requires all named parents to be exported samples in
#' the same family. Parent sex, self-parenting and ancestry cycles are checked.
#' At least one parent-offspring relationship is required. Unsequenced parents
#' and other more general Lep-MAP3 inputs are not supported by this exporter.
#' \item Pedigree rows are family, individual, father, mother, sex, phenotype;
#' each starts with CHR and POS. The command is a template, not an executed
#' analysis. No sex-chromosome model, segregation filter or linkage threshold
#' is selected automatically. Population HWE/LD filters should not be applied
#' indiscriminately to mapping families.
#' }
#' @references Rastas, P. (2017). Lep-MAP3: robust linkage mapping even for
#' low-coverage whole genome sequencing data. Bioinformatics, 33, 3726-3732.
#' \doi{10.1093/bioinformatics/btx494}.
#' @seealso \href{https://sourceforge.net/p/lep-map3/wiki/LM3%20Home/}{Lep-MAP3
#' input specification}, \code{read_vcf}, \code{write_vcf}.
#' @examples
#' \dontrun{
#' pedigree <- data.frame(FAMILY = "F1",
#'   INDIVIDUALS = c("father", "mother", "offspring"),
#'   FATHER = c("0", "0", "father"), MOTHER = c("0", "0", "mother"),
#'   SEX = c(1, 2, 0))
#' write_lepmap3("family.gds", pedigree, filename = "family")
#' }
#' @author Thierry Gosselin \email{thierrygosselin@@icloud.com}
#' @export
write_lepmap3 <- function(data, pedigree, likelihood = c("auto", "PL", "GL"),
  filename = NULL, path.folder = getwd(), overwrite = FALSE, verbose = TRUE) {
  likelihood <- match.arg(likelihood)
  .export_flag(overwrite, "overwrite")
  .export_flag(verbose, "verbose")
  if (verbose) .export_message("genometranslator::write_lepmap3")
  context <- .export_gds(data)
  on.exit(.export_close(context), add = TRUE)
  gds <- context$gds
  ids <- as.character(SeqArray::seqGetData(gds, "sample.id"))
  ped <- .lepmap3_pedigree(pedigree, ids)
  alleles <- SeqArray::seqGetData(gds, "allele")
  if (any(!grepl("^[ACGT],[ACGT]$", alleles)) ||
      any(substr(alleles, 1, 1) == substr(alleles, 3, 3)))
    stop("Lep-MAP3 export requires biallelic A/C/G/T SNPs.")
  fields <- c("PL", "GL")
  present <- vapply(fields, function(f)
    !is.null(gdsfmt::index.gdsn(gds, paste0("annotation/format/", f),
      silent = TRUE)), logical(1))
  if (likelihood == "auto") {
    if (!any(present)) stop("No PL or GL likelihoods retained in this GDS.")
    likelihood <- fields[which(present)[1L]]
  }
  if (!present[[likelihood]]) stop("Requested likelihood field is absent: ", likelihood)
  paths <- .export_paths(filename, "lepmap3", path.folder,
    c(".vcf", ".pedigree.tsv", ".samples.tsv", ".markers.tsv", ".command.txt"),
    overwrite)
  names(paths) <- c("vcf", "pedigree", "samples", "markers", "command")
  stage <- tempfile("lepmap3-")
  dir.create(stage)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)
  staged <- file.path(stage, basename(paths))
  raw <- tempfile("raw-", tmpdir = stage, fileext = ".vcf")
  SeqArray::seqGDS2VCF(gds, raw, verbose = FALSE)
  .lepmap3_vcf(raw, staged[1], ids, likelihood)
  utils::write.table(cbind("CHR", "POS", t(as.matrix(ped))), staged[2],
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  readr::write_tsv(tibble::as_tibble(ped), staged[3])
  mm <- genometranslator::extract_markers_metadata(gds, whitelist = FALSE)
  variants <- SeqArray::seqGetData(gds, "variant.id")
  mm <- mm[match(variants, mm$VARIANT_ID), , drop = FALSE]
  markers <- if ("MARKERS" %in% names(mm)) mm$MARKERS else
    rep(NA_character_, length(variants))
  readr::write_tsv(tibble::tibble(ROW = seq_along(variants),
    VARIANT_ID = variants, MARKERS = markers,
    VCF_ID = SeqArray::seqGetData(gds, "annotation/id"),
    CHROM = SeqArray::seqGetData(gds, "chromosome"),
    POS = SeqArray::seqGetData(gds, "position"),
    REF = substr(alleles, 1, 1), ALT = substr(alleles, 3, 3)), staged[4])
  writeLines(c("# Replace /path/to/Lep-MAP3/bin with your installation path.",
    "# Review ParentCall2 options for your pedigree and study design.",
    paste("java -cp /path/to/Lep-MAP3/bin ParentCall2",
      shQuote(paste0("data=", paths[2])),
      shQuote(paste0("vcfFile=", paths[1])),
      ">", shQuote(paste0(paths[1], ".parentcall.txt")))), staged[5])
  .export_publish(staged, paths, overwrite)
  if (verbose) .export_message("Exported ", length(ids), " samples and ",
    length(variants), " SNPs using ", likelihood, ".\nOutput: ", paths[1])
  invisible(paths)
}

.lepmap3_pedigree <- function(pedigree, ids) {
  if (is.character(pedigree) && length(pedigree) == 1L)
    pedigree <- readr::read_tsv(pedigree,
      col_types = readr::cols(.default = "c"), progress = FALSE)
  columns <- c("FAMILY", "INDIVIDUALS", "FATHER", "MOTHER", "SEX")
  if (!is.data.frame(pedigree) || !all(columns %in% names(pedigree)))
    stop("Pedigree requires FAMILY, INDIVIDUALS, FATHER, MOTHER and SEX.")
  if (!"PHENOTYPE" %in% names(pedigree)) pedigree$PHENOTYPE <- "0"
  p <- as.data.frame(lapply(pedigree[c(columns, "PHENOTYPE")], as.character),
    stringsAsFactors = FALSE)
  if (anyNA(p) || any(!nzchar(as.matrix(p))) ||
      any(grepl("[[:space:]#]", as.matrix(p))))
    stop("Pedigree values cannot be missing, empty, or contain whitespace/#.")
  if (anyDuplicated(p$INDIVIDUALS) || any(p$INDIVIDUALS == "0") ||
      !setequal(p$INDIVIDUALS, ids))
    stop("Pedigree IDs must be unique, nonzero and match active samples exactly.")
  p <- p[match(ids, p$INDIVIDUALS), , drop = FALSE]
  if (any(!p$SEX %in% c("0", "1", "2"))) stop("SEX must be 0, 1 or 2.")
  edges <- matrix(NA_integer_, nrow(p), 2)
  for (j in 1:2) {
    parent <- p[[c("FATHER", "MOTHER")[j]]]
    known <- parent != "0"
    ix <- match(parent[known], ids)
    if (anyNA(ix)) stop("Every named parent must be an active exported sample.")
    if (any(parent[known] == ids[known])) stop("A sample cannot be its own parent.")
    if (any(p$FAMILY[ix] != p$FAMILY[known])) stop("Parents must share the child's FAMILY.")
    if (any(p$SEX[ix] != as.character(j))) stop("Named parents require matching SEX (father 1, mother 2).")
    edges[known, j] <- ix
  }
  if (all(is.na(edges))) stop("Pedigree requires a parent-offspring relationship.")
  remaining <- seq_len(nrow(p))
  while (length(remaining)) {
    roots <- remaining[apply(edges[remaining, , drop = FALSE], 1,
      function(x) !any(x %in% remaining))]
    if (!length(roots)) stop("Pedigree contains an ancestry cycle.")
    remaining <- setdiff(remaining, roots)
  }
  p
}

# Stream the temporary VCF; never materialize the full genotype array.
.lepmap3_vcf <- function(source, destination, ids, likelihood) {
  input <- file(source, "rt"); on.exit(close(input), add = TRUE)
  output <- file(destination, "wt"); on.exit(close(output), add = TRUE)
  header <- FALSE
  repeat {
    lines <- readLines(input, n = 256L, warn = FALSE)
    if (!length(lines)) break
    for (line in lines) {
      if (startsWith(line, "##FORMAT=")) {
        if (grepl(paste0("^##FORMAT=<ID=(GT|", likelihood, "),"), line))
          writeLines(line, output)
      } else if (startsWith(line, "##")) {
        if (!startsWith(line, "##INFO=")) writeLines(line, output)
      } else if (startsWith(line, "#CHROM")) {
        if (!identical(strsplit(line, "\t", fixed = TRUE)[[1]][-(1:9)], ids))
          stop("Exported VCF sample order differs from the pedigree.")
        header <- TRUE
        writeLines(line, output)
      } else {
        f <- strsplit(line, "\t", fixed = TRUE)[[1]]
        if (!header || length(f) != 9L + length(ids)) stop("Malformed exported VCF.")
        fmt <- strsplit(f[9], ":", fixed = TRUE)[[1]]
        indexes <- match(c("GT", likelihood), fmt)
        if (anyNA(indexes)) stop("Every marker must contain GT and ", likelihood, ".")
        observed <- 0L
        for (i in seq_along(ids)) {
          values <- strsplit(f[9L+i], ":", fixed = TRUE)[[1]]
          gt <- values[indexes[1]]
          if (is.na(gt) || !grepl("^[01.][/|][01.]$", gt))
            stop("Diploid biallelic GT required at ", f[1], ":", f[2], ".")
          v <- values[indexes[2]]
          if (is.na(v) || v %in% c(".", ".,.,.")) v <- "." else {
            numbers <- suppressWarnings(as.numeric(strsplit(v, ",", fixed = TRUE)[[1]]))
            if (length(numbers) != 3L || any(!is.finite(numbers)) ||
                (likelihood == "PL" && any(numbers < 0 | numbers != floor(numbers) |
                  numbers > .Machine$integer.max)))
              stop("Invalid ", likelihood, " vector at ", f[1], ":", f[2],
                " for ", ids[i], ".")
            observed <- observed + 1L
          }
          f[9L+i] <- paste(gt, v, sep = ":")
        }
        if (!observed) stop("Marker has no observed ", likelihood, ": ", f[1], ":", f[2])
        f[8] <- "."; f[9] <- paste("GT", likelihood, sep = ":")
        writeLines(paste(f, collapse = "\t"), output)
      }
    }
  }
}
