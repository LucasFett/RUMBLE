#' Prepare Microbiome Data for Machine Learning
#'
#' Performs standard preprocessing for biomarker analysis: optional taxonomic
#' aggregation, optional removal of unclassified taxa, prevalence filtering,
#' NA imputation, and centered log-ratio (CLR) or robust CLR (rCLR) transformation.
#'
#' @param input A \code{phyloseq} object OR a numeric matrix/data.frame
#'   with samples as rows and taxa as columns.
#' @param metadata A data.frame containing sample metadata. Required when
#'   \code{input} is a matrix or data.frame. Must have rownames matching
#'   the sample identifiers in \code{input}.
#' @param tax_level Character. Taxonomic rank for aggregation (e.g.,
#'   \code{"Genus"}, \code{"Family"}). If \code{NULL} (default), no
#'   aggregation is performed.
#' @param min_prevalence Numeric. Minimum fraction of samples in which a
#'   taxon must be present to be retained (default 0.05, i.e., 5 percent).
#' @param min_abundance Numeric. Minimum relative abundance threshold to
#'   consider a taxon as present in a sample (default 0.0001).
#' @param pseudocount Numeric. Small value added before CLR transformation
#'   to handle zeros (default 1e-6).
#' @param remove_unclassified Logical. If \code{TRUE}, removes taxa whose
#'   names contain patterns such as \code{"uncultured"}, \code{"unknown"},
#'   \code{"unclassified"}, or \code{"_"}. Default is \code{FALSE}.
#' @param normalization_method Character. Method for compositional normalization.
#'   Options are \code{"clr"} (Centered Log-Ratio, default) or \code{"rclr"}
#'   (Robust CLR). The rCLR method preserves the zero structure of the data
#'   and is more robust to pseudocount artifacts, making it suitable for
#'   tree-based algorithms that may overfit on pseudocount noise. Default is
#'   \code{"clr"} to maintain compatibility with classical microbiome literature.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages.
#'
#' @return A named list with four elements:
#' \itemize{
#'   \item \code{features}: A data.frame of normalized abundances (samples x taxa).
#'   \item \code{metadata}: A data.frame of sample metadata aligned with features.
#'   \item \code{n_taxa_removed}: Integer. Total number of taxa removed.
#'   \item \code{n_samples}: Integer. Number of samples in the final dataset.
#' }
#'
#' @details
#' The preprocessing pipeline follows these steps:
#' \enumerate{
#'   \item Convert input to a \code{phyloseq} object (if not already).
#'   \item Aggregate taxa at the specified taxonomic level (optional).
#'   \item Remove unclassified/invalid taxa names (optional).
#'   \item Sanitise taxa names with \code{make.names()} for ML compatibility.
#'   \item Impute \code{NA} values with zero.
#'   \item Convert to relative abundance and filter by prevalence.
#'   \item Apply CLR or rCLR transformation.
#' }
#'
#' @examples
#'   # Load example data provided by the package
#'   ps_path <- system.file("extdata", "PRJEB38465_phyloseq_com_metadados_completos.rds",
#'                          package = "RUMBLE")
#'
#'   # Check if file exists (just in case)
#'   if (file.exists(ps_path)) {
#'     ps <- readRDS(ps_path)
#'
#'     # Run preprocessing with default CLR
#'     prep <- prepareData(ps, tax_level = "Genus", min_prevalence = 0.1)
#'
#'     # Or use robust CLR for tree-based models
#'     prep_rclr <- prepareData(ps, tax_level = "Genus",
#'                              normalization_method = "rclr")
#'
#'     # Check results
#'     head(prep$features[, 1:5])
#'   }
#'
#'
#' @importFrom methods is
#' @importFrom phyloseq phyloseq otu_table sample_data tax_table
#'      taxa_are_rows taxa_names rank_names prune_taxa ntaxa nsamples
#' @importFrom microbiome aggregate_taxa transform
#' @export
prepareData <- function(input,
                        metadata = NULL,
                        tax_level = NULL,
                        min_prevalence = 0.05,
                        min_abundance = 0.0001,
                        pseudocount = 1e-6,
                        remove_unclassified = FALSE,
                        normalization_method = "clr",
                        verbose = TRUE) {

  msg <- function(...) {
    if (verbose) message(...)
  }

  ## ------------------------------------------------------------------
  ## Validate normalization_method
  ## ------------------------------------------------------------------
  normalization_method <- match.arg(normalization_method,
                                    c("clr", "rclr"))
  msg("Normalization method: ", normalization_method)

  ## ------------------------------------------------------------------
  ## 1. Normalise input to phyloseq
  ## ------------------------------------------------------------------
  if (methods::is(input, "phyloseq")) {
    ps <- input
    if (is.null(metadata)) {
      metadata <- as.data.frame(phyloseq::sample_data(ps))
    }
  } else if (is.matrix(input) || is.data.frame(input)) {
    if (is.null(metadata)) {
      stop("'metadata' is required when 'input' is a ",
           "matrix or data.frame.")
    }
    input <- as.data.frame(input)
    # Simple NA imputation for input matrix
    input[is.na(input)] <- 0

    samples_common <- intersect(rownames(input), rownames(metadata))
    if (length(samples_common) == 0L) {
      stop("No matching sample names (rownames) between ",
           "'input' and 'metadata'.")
    }

    input_filt <- input[samples_common, , drop = FALSE]
    meta_filt  <- metadata[samples_common, , drop = FALSE]

    ps <- phyloseq::phyloseq(
      phyloseq::otu_table(as.matrix(input_filt),
                          taxa_are_rows = FALSE),
      phyloseq::sample_data(meta_filt)
    )
  } else {
    stop("Invalid input. Provide a 'phyloseq' object or a ",
         "matrix/data.frame with samples as rows.")
  }

  n_taxa_initial <- phyloseq::ntaxa(ps)
  msg("Input: ", phyloseq::nsamples(ps), " samples, ",
      n_taxa_initial, " taxa")

  ## ------------------------------------------------------------------
  ## 2. Taxonomic aggregation (optional)
  ## ------------------------------------------------------------------
  if (!is.null(tax_level)) {
    tt <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
    if (is.null(tt)) {
      warning("Aggregation requested but no tax_table found. ",
              "Skipping aggregation.")
    } else {
      ranks <- phyloseq::rank_names(ps)
      if (!tax_level %in% ranks) {
        stop("Taxonomic level '", tax_level,
             "' not found. Available: ",
             paste(ranks, collapse = ", "))
      }
      msg("Aggregating to: ", tax_level)
      ps <- microbiome::aggregate_taxa(ps, level = tax_level)
      msg("Taxa after aggregation: ", phyloseq::ntaxa(ps))
    }
  }

  ## ------------------------------------------------------------------
  ## 3. Remove unclassified / invalid taxa (optional)
  ## ------------------------------------------------------------------
  if (remove_unclassified) {
    taxa_nomes <- phyloseq::taxa_names(ps)
    invalid_pattern <- "uncultured|unknown|unclassified|_"
    valid_taxa <- taxa_nomes[
      !grepl(invalid_pattern, taxa_nomes, ignore.case = TRUE)
    ]
    valid_taxa <- valid_taxa[!is.na(valid_taxa) & valid_taxa != ""]

    n_removed_unclass <- length(taxa_nomes) - length(valid_taxa)

    if (n_removed_unclass > 0L) {
      msg("Removed ", n_removed_unclass,
          " unclassified/invalid taxa")
      if (length(valid_taxa) > 0) {
        ps <- phyloseq::prune_taxa(valid_taxa, ps)
      } else {
        stop("All taxa were filtered out by 'remove_unclassified'.")
      }
    }
  }

  ## ------------------------------------------------------------------
  ## 4. Sanitise taxa names for ML compatibility
  ## ------------------------------------------------------------------
  phyloseq::taxa_names(ps) <- make.names(
    phyloseq::taxa_names(ps), unique = TRUE
  )

  ## ------------------------------------------------------------------
  ## 5. Impute NAs with zero (on phyloseq object)
  ## ------------------------------------------------------------------
  otu_mat <- as.matrix(phyloseq::otu_table(ps))
  if (any(is.na(otu_mat))) {
    otu_mat[is.na(otu_mat)] <- 0
    phyloseq::otu_table(ps) <- phyloseq::otu_table(
      otu_mat, taxa_are_rows = phyloseq::taxa_are_rows(ps)
    )
  }

  ## ------------------------------------------------------------------
  ## 6. Prevalence filtering
  ## ------------------------------------------------------------------
  # Convert to compositional just for prevalence calculation
  ps_rel <- microbiome::transform(ps, "compositional")
  otu_rel <- as.matrix(phyloseq::otu_table(ps_rel))
  if (phyloseq::taxa_are_rows(ps_rel)) otu_rel <- t(otu_rel)

  # Check abundance against threshold
  prevalence <- vapply(
    seq_len(ncol(otu_rel)),
    function(j) mean(otu_rel[, j] > min_abundance),
    numeric(1L)
  )
  names(prevalence) <- colnames(otu_rel)
  keep_taxa <- names(prevalence)[prevalence >= min_prevalence]

  n_removed_prev <- phyloseq::ntaxa(ps) - length(keep_taxa)
  msg("Prevalence filter (>= ", min_prevalence * 100,
      "%): removed ", n_removed_prev, " taxa")

  if (length(keep_taxa) < 2L) {
    stop("Filtering too strict: fewer than 2 taxa remaining. ",
         "Reduce 'min_prevalence' or 'min_abundance'.")
  }
  ps_filtered <- phyloseq::prune_taxa(keep_taxa, ps)

  ## ------------------------------------------------------------------
  ## 7. Compositional normalization (CLR or rCLR)
  ## ------------------------------------------------------------------
  if (normalization_method == "clr") {
    # Standard CLR transformation
    ps_clr <- microbiome::transform(
      ps_filtered, "clr", pseudocount = pseudocount
    )
    msg("Applied CLR transformation with pseudocount = ", pseudocount)

  } else if (normalization_method == "rclr") {
    # Robust CLR (rCLR): Preserves zero structure
    otu_mat <- as(phyloseq::otu_table(ps_filtered), "matrix")
    taxa_rows <- phyloseq::taxa_are_rows(ps_filtered)

    if (taxa_rows) otu_mat <- t(otu_mat)

    # Aplica a transformação especificando MARGIN = 1 (amostras como linhas)
    otu_rclr <- vegan::decostand(otu_mat, method = "rclr", MARGIN = 1)

    # (Mata qualquer ruído de ponto flutuante do R)
    otu_rclr[otu_mat == 0] <- 0

    rownames(otu_rclr) <- rownames(otu_mat)
    colnames(otu_rclr) <- colnames(otu_mat)

    if (taxa_rows) otu_rclr <- t(otu_rclr)

    phyloseq::otu_table(ps_filtered) <- phyloseq::otu_table(otu_rclr, taxa_are_rows = taxa_rows)
    ps_clr <- ps_filtered

    msg("Applied rCLR (Robust CLR) transformation preserving zero structure (via vegan)")
  }

  ## ------------------------------------------------------------------
  ## 8. Build output
  ## ------------------------------------------------------------------
  features <- as.data.frame(phyloseq::otu_table(ps_clr))
  if (phyloseq::taxa_are_rows(ps_clr)) {
    features <- as.data.frame(t(features))
  }
  final_metadata <- as.data.frame(phyloseq::sample_data(ps_clr))

  n_taxa_final <- ncol(features)
  n_taxa_total_removed <- n_taxa_initial - n_taxa_final

  msg("Done: ", nrow(features), " samples, ",
      n_taxa_final, " features (",
      n_taxa_total_removed, " taxa removed)")

  list(
    features = features,
    metadata = final_metadata,
    n_taxa_removed = n_taxa_total_removed,
    n_samples = nrow(features)
  )
}
