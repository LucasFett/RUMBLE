#' Run the RUMBLE Pipeline
#'
#' Executes the complete RUMBLE analysis pipeline with publication-quality
#' settings. This function performs data preparation, model training with
#' hyperparameter tuning, evaluation on a held-out test set, and model-agnostic
#' interpretability via SHAP values and permutation importance.
#'
#' \strong{Default Parameters:} The default values are configured for
#' publication-quality analysis. All parameters are fully adjustable to allow
#' for exploratory analysis or more conservative/aggressive tuning as needed.
#'
#' All results, including plots and tables, are returned in a structured list
#' and optionally saved to disk.
#'
#' @param input A \code{phyloseq} object OR a numeric matrix/data.frame
#'   with samples as rows and taxa as columns.
#' @param metadata A data.frame of sample metadata. Required when
#'   \code{input} is a matrix or data.frame. Must have rownames
#'   matching the sample identifiers in \code{input}.
#' @param outcome_var Character. Name of the binary outcome variable in
#'   the metadata (e.g., \code{"ses"}).
#' @param class_of_interest Character. Label of the class of interest for
#'   disease/condition association analysis (e.g., \code{"CRC"}, \code{"disease"}).
#'   This parameter is mandatory and must match one of the levels in the
#'   \code{outcome_var} variable. SHAP values and feature directions will be
#'   calculated with respect to this class. For example, if
#'   \code{class_of_interest = "CRC"}, positive SHAP values indicate features
#'   that increase the probability of CRC, while negative values indicate
#'   protective features.
#' @param tax_level Character or \code{NULL}. Taxonomic rank for
#'   aggregation (e.g., \code{"Genus"}). Set to \code{NULL} to skip
#'   aggregation (default).
#' @param project_name Character. A label for the analysis, used in
#'   output file names (default \code{"RUMBLE_analysis"}).
#' @param output_dir Character. Directory path where output files
#'   (plots, tables) are saved. If \code{NULL}, no files are written
#'   (default \code{NULL}).
#' @param n_cores Integer. Number of CPU cores for parallel
#'   hyperparameter tuning (default 1).
#' @param train_prop Numeric. Proportion of data used for training
#'   (default 0.7).
#' @param cv_folds Integer. Number of cross-validation folds for tuning
#'   (default 5).
#' @param grid_size Integer. Number of hyperparameter combinations to
#'   evaluate per model during hyperparameter tuning (default 30).
#'   Higher values increase computation time but provide more stable
#'   hyperparameter selection. The default value (30) is recommended for
#'   publication-quality results. For quick exploratory analysis, you may
#'   reduce this to 10. For Meta-Analysis or Critical Study, consider 50 or higher.
#' @param top_n Integer. Number of top features to display in
#'   interpretability plots (default 20).
#' @param shap_reps Integer. Number of repetitions for SHAP and
#'   permutation importance calculations (default 100).
#'   Higher values provide more stable feature importance rankings.
#'   The default value (100) is recommended for publication-quality results.
#'   Feature rankings may vary by ±5% between runs at this level.
#'   For quick exploratory analysis, you may reduce this to 10 (feature rankings
#'   may vary by ±20%). For critical studies or meta-analyses, consider 1000
#'   repetitions (feature rankings stable to ±1%).
#'   \strong{Note:} This parameter significantly affects computation time.
#' @param xgb_trees Integer. Number of trees for the XGBoost model
#'   (default 1000). XGBoost is sensitive to the number of boosting iterations.
#'   The default value (1000) is recommended for publication-quality results.
#'   For quick exploratory analysis, you may reduce this to 500 to speed up
#'   computation. Beyond 1000, improvements are marginal and risk of overfitting
#'   increases.
#'   \strong{Note:} This parameter only affects the XGBoost model, not Random
#'   Forest, Elastic Net, or KNN models.
#' @param rf_trees Integer. Number of trees for the Random Forest model
#'   (default 500). The default value is well-established in the literature
#'   and appropriate for both exploratory and publication-quality analysis.
#'   Typically does not need adjustment.
#'   \strong{Note:} This parameter only affects the Random Forest model, not
#'   XGBoost, Elastic Net, or KNN models.
#' @param min_prevalence Numeric. Minimum prevalence threshold for
#'   taxa filtering (default 0.05).
#' @param min_abundance Numeric. Minimum relative abundance to
#'   consider a taxon present (default 0.0001).
#' @param remove_unclassified Logical. Whether to remove taxa with
#'   names matching common unclassified patterns (default
#'   \code{FALSE}).
#' @param seed Integer. Random seed for reproducibility (default 42).
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#'
#' @return A named list with five elements:
#'   \itemize{
#'     \item \code{models}: Named list of fitted model objects and predictions.
#'     \item \code{metrics}: Data.frame summarising performance metrics.
#'     \item \code{importance}: List containing SHAP and permutation results.
#'     \item \code{plots}: List of ggplot objects (ROC, Confusion Matrix, SHAP).
#'     \item \code{data}: List with train and test data.frames.
#'   }
#'
#' @examples
#'   # Load example data provided by the package
#'   ps_path <- system.file("extdata",
#'                          "PRJEB38465_phyloseq_com_metadados_completos.rds",
#'                          package = "RUMBLE")
#'
#'   if (file.exists(ps_path)) {
#'     ps <- readRDS(ps_path)
#'
#'     # EXAMPLE 1: Quick exploratory analysis Run pipeline with minimal parameters for demonstration speed
#'     # In real analysis, use defaults
#'     results <- RUMBLE(
#'       input = ps,
#'       outcome_var = "ses",
#'       class_of_interest = "Low",
#'       tax_level = "Genus",
#'       project_name = "Demo_Analysis",
#'       output_dir = tempdir(),
#'       cv_folds = 2,    # Reduced for example speed
#'       grid_size = 2,   # Reduced for example speed
#'       shap_reps = 1,   # Reduced for example speed
#'       n_cores = 1,
#'       verbose = FALSE
#'     )
#'
#'     # Inspect results
#'     print(results$metrics)
#'     print(results$plots$roc)
#'   }
#'     # EXAMPLE 2: Publication-quality analysis (default parameters)
#'     # Default parameters are configured for robust, publication-ready results
#'     results_publication <- RUMBLE(
#'       input = ps,
#'       outcome_var = "ses",
#'       class_of_interest = "Low",
#'       tax_level = "Genus",
#'       project_name = "Publication_Analysis",
#'       output_dir = tempdir(),
#'       # Using default parameters:
#'       # grid_size = 30 (robust hyperparameter tuning)
#'       # shap_reps = 100 (stable feature importance)
#'       # xgb_trees = 1000 (optimal XGBoost performance)
#'       n_cores = 1,
#'       verbose = TRUE
#'     )
#'
#'     # Inspect results
#'     print(results_publication$metrics)
#'     print(results_publication$plots$roc)
#'
#'
#' @importFrom readr write_tsv
#' @importFrom ggplot2 ggsave
#' @export
RUMBLE <- function(input,
                   metadata = NULL,
                   outcome_var,
                   class_of_interest,
                   tax_level = NULL,
                   project_name = "RUMBLE_analysis",
                   output_dir = NULL,
                   n_cores = 1L,
                   train_prop = 0.7,
                   cv_folds = 5L,
                   grid_size = 30L,
                   top_n = 20L,
                   shap_reps = 100L,
                   xgb_trees = 1000L,
                   rf_trees = 500L,
                   min_prevalence = 0.05,
                   min_abundance = 0.0001,
                   remove_unclassified = FALSE,
                   seed = 42L,
                   verbose = TRUE,
                   shap_data = "test") {

  msg <- function(...) {
    if (verbose) message(...)
  }

  ## ==================================================================
  ## 1. Data preparation
  ## ==================================================================
  msg("======================================================")
  msg("RUMBLE Pipeline: ", project_name)
  msg("======================================================")
  msg("[1/5] Preparing data...")

  prep <- prepareData(
    input = input,
    metadata = metadata,
    tax_level = tax_level,
    min_prevalence = min_prevalence,
    min_abundance = min_abundance,
    remove_unclassified = remove_unclassified,
    verbose = verbose
  )

  ## Build analysis data.frame
  if (!outcome_var %in% colnames(prep$metadata)) {
    stop("Variable '", outcome_var,
         "' not found in metadata. Available: ",
         paste(colnames(prep$metadata), collapse = ", "))
  }

  # Add outcome variable as a new column
  analysis_df <- prep$features
  analysis_df[[outcome_var]] <- prep$metadata[[outcome_var]]

  ## Remove samples with NA in outcome
  n_before <- nrow(analysis_df)
  analysis_df <- analysis_df[!is.na(analysis_df[[outcome_var]]), ,
                             drop = FALSE]
  n_after <- nrow(analysis_df)

  if (n_before > n_after) {
    msg("Removed ", n_before - n_after,
        " samples with NA in '", outcome_var, "'")
  }

  analysis_df[[outcome_var]] <- droplevels(
    as.factor(analysis_df[[outcome_var]])
  )

  lvls <- levels(analysis_df[[outcome_var]])
  if (length(lvls) != 2L) {
    stop("RUMBLE currently supports binary classification only. ",
         "Found ", length(lvls), " levels in '", outcome_var, "'.")
  }
  # Validate class_of_interest parameter
  if (missing(class_of_interest)) {
    stop("'class_of_interest' is required. Specify which class you want to study ",
         "(e.g., class_of_interest = 'CRC').")
  }

  if (!class_of_interest %in% lvls) {
    stop("'", class_of_interest, "' not found in '", outcome_var, "'. ",
         "Available classes: ", paste(lvls, collapse = ", "))
  }

  # Reorder factor levels to ensure class_of_interest is the second level
  # This ensures that SHAP values are calculated with respect to the class of interest
  negative_class <- lvls[lvls != class_of_interest]

  if (length(negative_class) != 1L) {
    stop("Expected exactly one negative class, but found ", length(negative_class), ".")
  }

  analysis_df[[outcome_var]] <- factor(
    analysis_df[[outcome_var]],
    levels = c(negative_class, class_of_interest)
  )

  msg("Reordered factor levels: '", negative_class, "' (reference) vs '",
      class_of_interest, "' (class of interest)")

  msg("Target: '", outcome_var, "' (",
      negative_class, " vs ", class_of_interest, ")")
  msg("Samples for analysis: ", nrow(analysis_df))

  ## ==================================================================
  ## 2. Train / test split
  ## ==================================================================
  msg("[2/5] Splitting data...")
  split <- .splitData(analysis_df, outcome_var,
                      prop = train_prop, seed = seed)
  train_data <- split$train
  test_data  <- split$test

  msg("Train: ", nrow(train_data), " | Test: ", nrow(test_data))

  ## ==================================================================
  ## 3. Model training & tuning
  ## ==================================================================
  msg("[3/5] Building and tuning models...")
  recipe <- .buildRecipe(train_data, outcome_var,
                         balance_classes = TRUE)
  workflows <- .buildWorkflows(recipe, xgb_trees = xgb_trees, rf_trees = rf_trees)

  tuned <- .tuneModels(
    workflows, train_data, outcome_var,
    grid_size = grid_size, cv_folds = cv_folds,
    n_cores = n_cores, seed = seed
  )

  ## ==================================================================
  ## 4. Evaluation
  ## ==================================================================
  msg("[4/5] Evaluating models on test set...")
  final_models <- .evaluateModels(
    tuned, train_data, test_data, outcome_var
  )

  ## ==================================================================
  ## 5. Interpretability
  ## ==================================================================
  msg("[5/5] Computing feature importance...")

  # Determine which data to use for SHAP
  if (shap_data == "test") {
    explainer_data <- test_data
    msg("Using TEST set for SHAP interpretation (default)")
  } else if (shap_data == "train") {
    explainer_data <- train_data
    msg("Using TRAIN set for SHAP interpretation")
  } else {
    stop("Invalid 'shap_data' parameter. Use 'test' or 'train'.")
  }

  explainers <- createExplainers(
    final_models, explainer_data, outcome_var,
    class_of_interest = class_of_interest
  )
  importance <- computeFeatureImportance(
    explainers, explainer_data, outcome_var,
    top_n = top_n, repetitions = shap_reps,
    n_cores = n_cores, verbose = verbose
  )

  ## ==================================================================
  ## Generate plots
  ## ==================================================================
  msg("Generating plots...")
  plots <- list(
    metrics   = plotMetricsComparison(final_models),
    roc       = plotRocCurves(final_models, outcome_var),
    cm        = plotConfusionMatrices(final_models, outcome_var),
    shap_spearman = plotShapGlobal(
      importance$global_importance$spearman, outcome_var,
      class_of_interest = class_of_interest,
      negative_class = negative_class,
      top_n = top_n,
      metric_name = "Spearman"
    ),
    shap_mean = plotShapGlobal(
      importance$global_importance$mean_shap, outcome_var,
      class_of_interest = class_of_interest,
      negative_class = negative_class,
      top_n = top_n,
      metric_name = "Mean SHAP"
    ),
    shap_beeswarm = plotShapBeeswarm(
      importance$shap_raw, top_n = top_n
    ),
    shap_prevalence = plotTaxaPrevalence(
      explainer_data, importance$global_importance$spearman,
      target_var = outcome_var, top_n = top_n
    ),
    biomarker_integrated_spearman = plotBiomarkerIntegrated(
      explainer_data, importance$global_importance$spearman,
      target_var = outcome_var, top_n = top_n,
      class_of_interest = class_of_interest,
      negative_class = negative_class,
      metric_name = "Spearman"
    ),
    biomarker_integrated_mean = plotBiomarkerIntegrated(
      explainer_data, importance$global_importance$mean_shap,
      target_var = outcome_var, top_n = top_n,
      class_of_interest = class_of_interest,
      negative_class = negative_class,
      metric_name = "Mean SHAP"
    ),
    perm_heat = plotPermutationHeatmap(
      importance$permutation_top
    )
  )

  ## Metrics summary table
  metrics_summary <- summariseMetrics(final_models)

  ## ==================================================================
  ## Save outputs (optional)
  ## ==================================================================
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    prefix <- file.path(output_dir, project_name)
    msg("Saving outputs to: ", output_dir)

    ## Plots
    plot_names <- c("metrics", "roc", "cm", "shap_spearman", "shap_mean",
                    "shap_beeswarm", "shap_prevalence",
                    "biomarker_integrated_spearman", "biomarker_integrated_mean", "perm_heat")

    # Define sizes (width, height)
    plot_sizes <- list(
      metrics                       = c(8, 5),
      roc                           = c(8, 5),
      cm                            = c(10, 8),
      shap_spearman                 = c(9, 8),
      shap_mean                     = c(9, 8),
      shap_beeswarm                 = c(10, 8),
      shap_prevalence               = c(9, 8),
      biomarker_integrated_spearman = c(12, 9),
      biomarker_integrated_mean     = c(12, 9),
      perm_heat                     = c(9, 8)
    )

    for (pname in plot_names) {
      if (!is.null(plots[[pname]])) {
        sz <- plot_sizes[[pname]]
        ggplot2::ggsave(
          filename = paste0(prefix, "_", pname, ".png"),
          plot = plots[[pname]],
          width = sz[1L], height = sz[2L], dpi = 300
        )
      }
    }

    ## Tables
    readr::write_tsv(
      metrics_summary,
      paste0(prefix, "_model_metrics.tsv")
    )
    readr::write_tsv(
      importance$global_importance$spearman,
      paste0(prefix, "_shap_global_spearman.tsv")
    )
    readr::write_tsv(
      importance$global_importance$mean_shap,
      paste0(prefix, "_shap_global_mean_shap.tsv")
    )
    readr::write_tsv(
      importance$permutation_top,
      paste0(prefix, "_permutation_top.tsv")
    )
    msg("All outputs saved.")
  }

  ## ==================================================================
  ## Return
  ## ==================================================================
  msg("======================================================")
  msg("RUMBLE Pipeline complete.")
  msg("======================================================")

  list(
    models     = final_models,
    metrics    = metrics_summary,
    importance = importance,
    plots      = plots,
    data       = list(train = train_data, test = test_data)
  )
}
