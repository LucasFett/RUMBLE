#' Run Differential Abundance Analysis
#'
#' Performs differential abundance (DA) analysis using ANCOM-BC, ALDEx2, or Wilcoxon
#' methods. This function provides a unified wrapper around state-of-the-art
#' compositional data analysis methods, each designed to handle the inherent
#' compositionality of microbiome data while controlling for false discovery rate (FDR).
#'
#' @details
#' \strong{Wilcoxon (Fallback):} Runs a basic Wilcoxon Rank-Sum test on relative abundances.
#' Does not require Bioconductor dependencies. Useful as a baseline or fallback method.
#'
#' \strong{ANCOM-BC:} Implements the ANCOM-BC2 method, which uses a
#' bias-correction framework specifically designed for compositional data. It provides
#' robust control of false discovery rate and high statistical power for detecting
#' differential abundance signals. ANCOM-BC is recommended for most analyses.
#'
#' \strong{ALDEx2 (Optional):} Uses additive log-ratio (ALR) transformation with
#' Monte Carlo sampling and Dirichlet-multinomial modeling. ALDEx2 is more conservative
#' and is recommended when ultra-conservative results are needed or for validation purposes.
#'
#' All methods accept raw count data (not CLR-transformed) and handle compositionality
#' internally through their respective statistical frameworks.
#'
#' @param counts A numeric matrix or data.frame with samples as rows and taxa as columns.
#'   Should contain raw counts (not CLR-transformed). Rownames should match sample identifiers
#'   in \code{metadata}.
#' @param metadata A data.frame with sample metadata. Must have rownames matching
#'   the sample identifiers in \code{counts}.
#' @param target_var Character. Name of the binary outcome variable in \code{metadata}.
#' @param method Character. Differential abundance method to use. Options are
#'   \code{"wilcoxon"} (default), \code{"ancombc"}, or \code{"aldex2"}.
#' @param verbose Logical. If \code{TRUE} (default), prints progress messages.
#'
#' @return A data.frame with standardized columns:
#'   \itemize{
#'     \item \code{Taxon}: Character. Taxon/feature name.
#'     \item \code{logFC}: Numeric. Log2 fold change.
#'     \item \code{p_val}: Numeric. Raw p-value.
#'     \item \code{adj_p_val}: Numeric. Adjusted p-value (FDR-corrected).
#'   }
#'
#' @examples
#'   # Create example count data
#'   set.seed(123)
#'   counts <- matrix(rpois(100, lambda = 5), nrow = 10, ncol = 10)
#'   rownames(counts) <- paste0("Sample_", 1:10)
#'   colnames(counts) <- paste0("Taxon_", 1:10)
#'
#'   # Create metadata
#'   metadata <- data.frame(
#'     row.names = paste0("Sample_", 1:10),
#'     group = c(rep("Control", 5), rep("Disease", 5))
#'   )
#'
#'   # Run differential abundance with Wilcoxon (fallback method)
#'   da_results <- runDifferentialAbundance(
#'     counts = counts,
#'     metadata = metadata,
#'     target_var = "group",
#'     method = "wilcoxon"
#'   )
#'
#'   head(da_results)
#'
#' @importFrom stats p.adjust as.formula
#' @export
runDifferentialAbundance <- function(counts,
                                     metadata,
                                     target_var,
                                     method = c("wilcoxon", "ancombc", "aldex2"),
                                     verbose = TRUE) {

  method <- match.arg(method)

  msg <- function(...) {
    if (verbose) message(...)
  }

  msg("Running Differential Abundance Analysis (", method, ")...")

  # Validate inputs
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("'counts' must be a matrix or data.frame.")
  }

  if (nrow(metadata) != nrow(counts)) {
    stop("Number of rows in 'counts' and 'metadata' must match.")
  }

  if (!all(rownames(counts) == rownames(metadata))) {
    stop("Rownames of 'counts' and 'metadata' must match exactly.")
  }

  if (!target_var %in% colnames(metadata)) {
    stop("'target_var' (", target_var, ") not found in 'metadata'.")
  }

  counts <- as.matrix(counts)

  # Check for relative abundances
  if (any(counts > 0 & counts < 1) && method != "wilcoxon") {
    stop("ERROR: ANCOM-BC and ALDEx2 require raw counts (integers). Your data appear to be relative abundances. Please provide raw count data to use these methods, or switch to 'wilcoxon'.")
  }

  if (method %in% c("ancombc", "aldex2")) {
    counts <- round(counts)
    storage.mode(counts) <- "integer"
  }

  if (method == "ancombc") {
    result <- .run_ancombc(counts, metadata, target_var, verbose = verbose)
  } else if (method == "aldex2") {
    result <- .run_aldex2(counts, metadata, target_var, verbose = verbose)
  } else if (method == "wilcoxon") {
    result <- .run_wilcoxon(counts, metadata, target_var, verbose = verbose)
  }

  msg("Differential Abundance Analysis completed.")
  result
}


#' Internal: Run ANCOM-BC Analysis
#'
#' @keywords internal
.run_ancombc <- function(counts, metadata, target_var, verbose = TRUE) {

  msg <- function(...) {
    if (verbose) message(...)
  }

  # Check if required packages are available
  if (!requireNamespace("TreeSummarizedExperiment", quietly = TRUE)) {
    stop("Package 'TreeSummarizedExperiment' is required for ANCOM-BC analysis. ",
         "Install it with: BiocManager::install('TreeSummarizedExperiment')")
  }

  if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
    stop("Package 'ANCOMBC' is required for ANCOM-BC analysis. ",
         "Install it with: BiocManager::install('ANCOMBC')")
  }

  msg("  - Preparing TreeSummarizedExperiment object...")

  # Create TreeSummarizedExperiment object
  tse <- TreeSummarizedExperiment::TreeSummarizedExperiment(
    assays = list(counts = t(counts)),  # Features as rows, samples as columns
    colData = metadata
  )

  msg("  - Running ANCOM-BC2...")

  # Run ANCOM-BC2
  ancombc_result <- ANCOMBC::ancombc2(
    data = tse,
    assay_name = "counts",
    tax_level = NULL,  # Use original features
    fix_formula = target_var,
    rand_formula = NULL,
    p_adj_method = "BH",  # Benjamini-Hochberg FDR correction
    prv_cut = 0,  # No prevalence filtering (already done in prepareData)
    lib_cut = 0,  # No library size filtering
    s0_perc = 0.05,
    neg_lb = TRUE,
    verbose = FALSE
  )

  msg("  - Extracting and standardizing results...")

  # Extract results
  res_df <- ancombc_result$res

  # Correctly get the levels to avoid alphabetical ordering mismatch
  if(is.factor(metadata[[target_var]])) {
    target_levels <- levels(metadata[[target_var]])
  } else {
    target_levels <- unique(metadata[[target_var]])
  }

  if (length(target_levels) != 2) {
    stop("Target variable must be binary.")
  }

  # =================================================================
  # FIX
  # =================================================================
  base_coef <- paste0(target_var, target_levels[2])
  coef_col <- paste0("lfc_", base_coef)
  pval_col <- paste0("p_", base_coef)
  qval_col <- paste0("q_", base_coef)

  # Check if columns exist
  if (!coef_col %in% colnames(res_df)) {
    # Try alternative naming, but STRICTLY EXCLUDING the Intercept
    lfc_cols <- colnames(res_df)[grep("lfc", colnames(res_df), ignore.case = TRUE)]
    lfc_cols <- lfc_cols[!grepl("Intercept", lfc_cols, ignore.case = TRUE)]
    coef_col <- lfc_cols[1]

    p_cols <- colnames(res_df)[grep("^p_", colnames(res_df), ignore.case = TRUE)]
    p_cols <- p_cols[!grepl("Intercept", p_cols, ignore.case = TRUE)]
    pval_col <- p_cols[1]

    q_cols <- colnames(res_df)[grep("^q_", colnames(res_df), ignore.case = TRUE)]
    q_cols <- q_cols[!grepl("Intercept", q_cols, ignore.case = TRUE)]
    qval_col <- q_cols[1]

    if (is.na(coef_col) || is.na(pval_col) || is.na(qval_col)) {
      stop("Could not identify coefficient, p-value, or q-value columns in ANCOM-BC output.")
    }
  }

  # Standardize output
  result <- data.frame(
    Taxon = res_df$taxon,
    logFC = res_df[[coef_col]],
    p_val = res_df[[pval_col]],
    adj_p_val = res_df[[qval_col]],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  result
}


#' Internal: Run ALDEx2 Analysis
#'
#' @keywords internal
.run_aldex2 <- function(counts, metadata, target_var, verbose = TRUE) {

  msg <- function(...) {
    if (verbose) message(...)
  }

  # Check if required packages are available
  if (!requireNamespace("ALDEx2", quietly = TRUE)) {
    stop("Package 'ALDEx2' is required for ALDEx2 analysis. ",
         "Install it with: BiocManager::install('ALDEx2')")
  }

  msg("  - Preparing data for ALDEx2...")

  # Ensure target is evaluated as factor to get strict 0 and 1 mapping
  target_levels <- levels(factor(metadata[[target_var]]))
  if (length(target_levels) != 2) {
    stop("Target variable must be binary.")
  }

  group_vector <- factor(metadata[[target_var]])
  group_numeric <- as.numeric(group_vector) - 1  # 0 and 1

  msg("  - Running ALDEx2 analysis...")

  # Run ALDEx2 analysis
  aldex_result <- ALDEx2::aldex(
    reads = t(counts),  # Features as rows, samples as columns
    conditions = group_numeric,
    mc.samples = 128,
    test = "t",
    effect = TRUE,
    verbose = FALSE
  )

  msg("  - Extracting and standardizing results...")

  # Extract results and standardize
  result <- data.frame(
    Taxon = rownames(aldex_result),
    logFC = aldex_result$effect,
    p_val = aldex_result$we.ep,
    adj_p_val = aldex_result$we.eBH,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  result
}


#' Internal: Run Wilcoxon Analysis
#'
#' @keywords internal
.run_wilcoxon <- function(counts, metadata, target_var, verbose = TRUE) {

  msg <- function(...) {
    if (verbose) message(...)
  }

  msg("  - Running basic Wilcoxon Rank-Sum test on relative abundances (Fallback method)...")

  target_levels <- levels(factor(metadata[[target_var]]))
  neg_class <- target_levels[1]
  pos_class <- target_levels[2]

  # Convert counts to relative abundance for Wilcoxon
  depths <- rowSums(counts)
  depths[depths == 0] <- 1 # Prevent division by zero
  rel_abund <- sweep(counts, 1, depths, "/")
  pseudo <- 1e-6 # small pseudocount for logFC

  res_list <- lapply(colnames(rel_abund), function(tax) {
    x <- rel_abund[, tax]

    # Skip if variance is 0
    if (length(unique(x)) < 2) {
      p_val <- 1.0
    } else {
      wt <- suppressWarnings(stats::wilcox.test(x ~ metadata[[target_var]]))
      p_val <- wt$p.value
    }

    mean_neg <- mean(x[metadata[[target_var]] == neg_class])
    mean_pos <- mean(x[metadata[[target_var]] == pos_class])
    logfc <- log2((mean_pos + pseudo) / (mean_neg + pseudo))

    data.frame(Taxon = tax, logFC = logfc, p_val = p_val, stringsAsFactors = FALSE)
  })

  result <- do.call(rbind, res_list)
  result$adj_p_val <- stats::p.adjust(result$p_val, method = "BH")

  return(result)
}
