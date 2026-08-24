#' Run the Complete RUMBLE Analysis Pipeline
#'
#' Executes the complete RUMBLE analysis pipeline with publication-quality
#' settings. This function performs data preparation, model training with
#' hyperparameter tuning via nested cross-validation, evaluation on held-out
#' test folds, and model-agnostic interpretability via SHAP values and
#' permutation importance.
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
#'   calculated with respect to this class.
#' @param tax_level Character or \code{NULL}. Taxonomic rank for
#'   aggregation (e.g., \code{"Genus"}). Set to \code{NULL} to skip
#'   aggregation (default).
#' @param correlation_threshold Numeric or \code{NULL}. Threshold for removing
#'   highly correlated features during preprocessing. Default is \code{NULL}
#'   (no filtering), which is highly recommended for microbiome data to preserve
#'   biological consortia/guilds.
#' @param project_name Character. A label for the analysis, used in
#'   output file names (default \code{"RUMBLE_analysis"}).
#' @param output_dir Character. Directory path where output files
#'   (plots, tables) are saved. If \code{NULL}, no files are written
#'   (default \code{NULL}).
#' @param n_cores Integer. Number of CPU cores for parallel
#'   hyperparameter tuning (default 1).
#' @param outer_folds Integer. Number of outer cross-validation folds for
#'   testing (default 5).
#' @param cv_folds Integer. Number of cross-validation folds for tuning
#'   (default 5).
#' @param grid_size Integer. Number of hyperparameter combinations to
#'   evaluate per model during hyperparameter tuning (default 30).
#' @param top_n Integer. Number of top features to display in
#'   interpretability plots (default 20).
#' @param shap_reps Integer. Number of repetitions for SHAP and
#'   permutation importance calculations (default 25).
#' @param shap_method Character. Method for SHAP calculation. Options are
#'   \code{"exact"} (default) or \code{"fast"}.
#' @param shap_model Character vector controlling which model-specific SHAP
#'   profiles are generated. Use \code{"all"} (default) to generate for every
#'   model, \code{"consensus"} to skip model-specific plots, or provide
#'   specific model names (e.g., \code{c("RF", "XGB")}).
#' @param xgb_trees Integer. Number of trees for XGBoost (default 1000).
#' @param rf_trees Integer. Number of trees for Random Forest (default 500).
#' @param min_prevalence Numeric. Minimum prevalence threshold for
#'   taxa filtering (default 0.05).
#' @param min_abundance Numeric. Minimum relative abundance to
#'   consider a taxon present (default 0.0001).
#' @param pseudocount Numeric. Pseudocount added to exact zeros before the
#'   Centered Log-Ratio (CLR) transformation (default 1e-6). Re-exposed as a
#'   top-level parameter after having been hardcoded at the three
#'   \code{.apply_clr()} call sites inside the pipeline; see
#'   \code{vignette("RUMBLE_tutorial")} or the manuscript's Section 2.1.1 for
#'   the rationale behind the default value and its numerical effect on the
#'   separation between structural zeros and rare biological signal.
#' @param min_fold_frequency Numeric. Minimum fraction of outer folds in which
#'   a feature must appear in the top-N to be included in the consensus ranking
#'   (default 0.7, i.e., 70 percent of folds). This filters out unstable
#'   features that appear sporadically (only applied if \code{apply_stability_filter = TRUE}).
#' @param remove_unclassified Logical. Whether to remove taxa with
#'   names matching common unclassified patterns (default \code{FALSE}).
#' @param unclassified_patterns Character. Regex pattern for identifying
#'   unclassified taxa. Only used when \code{remove_unclassified = TRUE}.
#'   Default: \code{"uncultured|unknown|unclassified"}.

#' @param class_balance_method Character. Method for handling class imbalance during
#'   model training. Options are \code{"downsample"} (default, reduces majority
#'   class to match minority class size), \code{"class_weights"} (applies native
#'   per-class weighting at the engine level -- \code{ranger::class.weights} for
#'   RF and \code{xgboost::scale_pos_weight} for XGB -- preserving 100 percent of
#'   the training data; ENET and KNN have no native per-class weighting mechanism
#'   and are fit unweighted under this mode, see \code{.buildWorkflows()}), or
#'   \code{"none"} (no balancing, useful for already-balanced data).
#' @param apply_stability_filter Logical. Whether to filter features based on
#'   \code{min_fold_frequency} in the final biomarker table and plots (default \code{FALSE}).
#'   When \code{FALSE}, all features are retained; when \code{TRUE}, only features
#'   appearing in at least \code{min_fold_frequency} fraction of folds are included.
#' @param seed Integer. Random seed for reproducibility (default 42).
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @param metric_cutoffs Named numeric vector. Optional quality filter for
#'   models based on performance metrics. Default is \code{NULL} (no filtering).
#' @param run_da Logical. Whether to run differential abundance analysis
#'   (default \code{TRUE}).
#' @param da_method Character. Differential abundance method to use.
#'   Options are \code{"wilcoxon"} (default), \code{"ancombc"}, or
#'   \code{"aldex2"}.
#'
#' @return A named list with the following elements:
#'   \itemize{
#'     \item \code{models}: Named list of fitted model objects and predictions.
#'     \item \code{metrics}: Data.frame summarising performance metrics.
#'     \item \code{importance}: List containing SHAP and permutation results.
#'     \item \code{plots}: List of ggplot objects.
#'     \item \code{selected_models}: Character vector of models used for
#'       interpretability after applying \code{metric_cutoffs}.
#'     \item \code{da_results}: Data.frame with differential abundance results.
#'     \item \code{shap_overfitting}: Train vs test SHAP stability metrics.
#'     \item \code{model_consistency}: Inter-model agreement metrics.
#'     \item \code{shap_consistency}: Inter-fold SHAP stability metrics.
#'     \item \code{feature_frequency}: Feature selection frequency across folds.
#'     \item \code{hyperparameter_stability}: Hyperparameter variance summary.
#'   }
#'
#' @examples
#'   ps_path <- system.file("extdata",
#'                          "PRJEB38465_phyloseq_com_metadados_completos.rds",
#'                          package = "RUMBLE")
#'
#'   if (file.exists(ps_path)) {
#'     ps <- readRDS(ps_path)
#'
#'     results <- RUMBLE(
#'       input = ps,
#'       outcome_var = "ses",
#'       class_of_interest = "Low",
#'       tax_level = "Genus",
#'       project_name = "Demo_Analysis",
#'       output_dir = tempdir(),
#'       cv_folds = 2,
#'       grid_size = 2,
#'       shap_reps = 1,
#'       n_cores = 1,
#'       verbose = FALSE,
#'       outer_folds = 3,
#'       xgb_trees = 5,
#'       rf_trees = 5
#'     )
#'
#'     print(results$metrics)
#'   }
#'
#' @importFrom readr write_tsv
#' @importFrom ggplot2 ggsave
#' @export
RUMBLE <- function(input,
                   metadata = NULL,
                   outcome_var,
                   class_of_interest,
                   tax_level = NULL,
                   correlation_threshold = NULL,
                   project_name = "RUMBLE_analysis",
                   output_dir = NULL,
                   n_cores = 1L,
                   outer_folds = 5L,
                   cv_folds = 5L,
                   grid_size = 30L,
                   top_n = 20L,
                   shap_reps = 25L,
                   xgb_trees = 1000L,
                   rf_trees = 500L,
                   min_prevalence = 0.05,
                   min_abundance = 0.0001,
                   pseudocount = 1e-6,
                   min_fold_frequency = 0.7,
                   remove_unclassified = FALSE,
                   unclassified_patterns = "uncultured|unknown|unclassified",
                   class_balance_method = c("downsample", "class_weights", "none"),
                   apply_stability_filter = FALSE,
                   seed = 42L,
                   verbose = TRUE,
                   shap_method = "exact",
                   shap_model = "all",
                   metric_cutoffs = NULL,
                   run_da = TRUE,
                   da_method = c("wilcoxon", "ancombc", "aldex2")) {

  msg <- function(...) {
    if (verbose) message(...)
  }

  ## ==================================================================
  ## 0. EARLY VALIDATION & CLEANUP
  ## ==================================================================

  # Ensure parallel workers are cleaned up on exit, even in case of
  # error. This prevents zombie workers from lingering between RUMBLE()
  # calls, which can cause hangs.
  on.exit(future::plan("sequential"), add = TRUE)

  # Validate da_method parameter
  da_method <- match.arg(da_method)

  # Validate class_balance_method parameter
  class_balance_method <- match.arg(class_balance_method)

  # Early validation: Check if input is relative abundance when DA requires counts

  if (run_da && da_method != "wilcoxon") {
    if (is.matrix(input) || is.data.frame(input)) {
      input_check <- as.matrix(input)
      if (any(input_check > 0 & input_check < 1, na.rm = TRUE)) {
        stop("ANCOM-BC and ALDEx2 require raw integer counts, not relative abundances. ",
             "Your input contains values between 0 and 1. ",
             "Either provide raw counts or use da_method = 'wilcoxon'.")
      }
    }
  }

  # Validate shap_model parameter early
  valid_shap_model_keywords <- c("all", "consensus")
  if (!all(shap_model %in% valid_shap_model_keywords)) {
    # Will be validated after models are trained (against actual model names)
    shap_model_is_custom <- TRUE
  } else {
    shap_model_is_custom <- FALSE
  }

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

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
    min_prevalence = 0.0,
    min_abundance = min_abundance,
    remove_unclassified = remove_unclassified,
    unclassified_patterns = unclassified_patterns,
    verbose = verbose
  )

  # Store filtered_counts for later use in plotting and DA analysis
  filtered_counts <- prep$filtered_counts

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
  valid_idx <- !is.na(analysis_df[[outcome_var]])
  n_before <- nrow(analysis_df)

  analysis_df <- analysis_df[valid_idx, , drop = FALSE]
  filtered_counts <- filtered_counts[valid_idx, , drop = FALSE]
  prep$metadata <- prep$metadata[valid_idx, , drop = FALSE]

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
  negative_class <- lvls[lvls != class_of_interest]

  if (length(negative_class) != 1L) {
    stop("Expected exactly one negative class, but found ", length(negative_class), ".")
  }

  analysis_df[[outcome_var]] <- factor(
    analysis_df[[outcome_var]],
    levels = c(negative_class, class_of_interest)
  )

  prep$metadata[[outcome_var]] <- analysis_df[[outcome_var]]

  msg("Reordered factor levels: '", negative_class, "' (reference) vs '",
      class_of_interest, "' (class of interest)")
  msg("Target: '", outcome_var, "' (",
      negative_class, " vs ", class_of_interest, ")")
  msg("Samples for analysis: ", nrow(analysis_df))

  ## ==================================================================
  ## 1.5 Differential Abundance Analysis
  ## ==================================================================
  da_results <- NULL

  if (run_da) {
    msg("[1.5/5] Running differential abundance analysis (", da_method, ")...")

    da_results <- runDifferentialAbundance(
      counts = filtered_counts,
      metadata = prep$metadata,
      target_var = outcome_var,
      method = da_method,
      verbose = verbose
    )

    msg("Differential abundance analysis complete. Found ",
        nrow(da_results), " features.")
  } else {
    msg("Skipping differential abundance analysis (run_da = FALSE)")
  }

  ## ==================================================================
  ## 2. Nested Cross-Validation Setup
  ## ==================================================================
  msg("[2/5] Creating ", outer_folds, " outer folds for Nested CV...")
  outer_folds_obj <- .createOuterFolds(analysis_df, outcome_var, folds = outer_folds, seed = seed)

  ## ==================================================================
  ## 3 & 4. Train, Tune, and Evaluate across Folds
  ## ==================================================================
  msg("[3/5] & [4/5] Running Nested CV (Tuning and Evaluating)...")

  fold_results <- purrr::imap(outer_folds_obj$splits, function(split_obj, fold_idx) {
    fold_name <- paste0("Fold_", fold_idx)
    msg("\n=====================================")
    msg("  Processing ", fold_name)
    msg("=====================================")

    train_data <- rsample::training(split_obj)
    test_data  <- rsample::testing(split_obj)

    msg("  Train: ", nrow(train_data), " | Test: ", nrow(test_data))

    # ---------------------------------------------------------------
    #%% [ORIGINAL -- pre-commit1] CLR was applied BEFORE the intra-fold
    #%% prevalence filter, violating what Section 2.1.1 of the manuscript
    #%% promises. The per-sample geometric mean used in the CLR included
    #%% taxa that would be discarded immediately afterward by the filter.
    #%%
    #%% # ---------------------------------------------------------------
    #%% # STRICT SEPARATION: CLR applied intra-fold to avoid Data Leakage criticism
    #%% # ---------------------------------------------------------------
    #%% train_data <- .apply_clr(train_data, outcome_var, pseudocount = 1e-6)
    #%% test_data  <- .apply_clr(test_data, outcome_var, pseudocount = 1e-6)
    #%% # ---------------------------------------------------------------
    #%% # STRICT ANTI-LEAKAGE FILTER: Dynamic Prevalence based on Train only
    #%% # ---------------------------------------------------------------
    #%% train_ids <- rownames(train_data)
    #%% train_counts <- filtered_counts[rownames(filtered_counts) %in% train_ids, , drop = FALSE]
    #%% train_depths <- rowSums(train_counts)
    #%% train_depths[train_depths == 0] <- 1
    #%% train_rel <- sweep(train_counts, 1, train_depths, "/")
    #%% prev_train <- colMeans(train_rel > min_abundance)
    #%% valid_taxa <- names(prev_train)[prev_train >= min_prevalence]
    #%% cols_to_keep <- c(valid_taxa, outcome_var)
    #%% train_data_filtered <- train_data[, colnames(train_data) %in% cols_to_keep, drop = FALSE]
    #%% test_data_filtered  <- test_data[, colnames(test_data) %in% cols_to_keep, drop = FALSE]
    # ---------------------------------------------------------------
    # STRICT ANTI-LEAKAGE FILTER: Dynamic Prevalence based on Train only
    # Executed BEFORE the CLR (post-submission fix): Section 2.1.1 of the
    # manuscript promises that the CLR operates only on taxa that survive
    # the intra-fold prevalence filter. Before this fix, the CLR was
    # computed on the full (pre-filter) table, causing the per-sample
    # log-geometric-mean to include taxa that would be discarded next --
    # inconsistent with both the text and the anti-leakage design.
    # ---------------------------------------------------------------
    train_ids <- rownames(train_data)
    train_counts <- filtered_counts[rownames(filtered_counts) %in% train_ids, , drop = FALSE]

    # Calculate relative abundance internally for prevalence assessment
    train_depths <- rowSums(train_counts)
    train_depths[train_depths == 0] <- 1
    train_rel <- sweep(train_counts, 1, train_depths, "/")

    # Calculate prevalence looking only at Training data
    prev_train <- colMeans(train_rel > min_abundance)
    valid_taxa <- names(prev_train)[prev_train >= min_prevalence]

    # Keep only validated taxa + target variable (still pre-CLR here)
    cols_to_keep <- c(valid_taxa, outcome_var)

    train_data_filtered <- train_data[, colnames(train_data) %in% cols_to_keep, drop = FALSE]
    test_data_filtered  <- test_data[, colnames(test_data) %in% cols_to_keep, drop = FALSE]

    msg("  Anti-Leakage: Retained ", length(valid_taxa), " taxa prevalent in training fold.")
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    # STRICT SEPARATION: CLR applied intra-fold, now ONLY on the taxa that
    # survived the prevalence filter above. The training set defines the
    # retained taxa (train_data_filtered); the same set of columns is
    # applied to the test set (test_data_filtered) before the CLR,
    # ensuring both sets share the same compositional space.
    # ---------------------------------------------------------------
    train_data_filtered <- .apply_clr(train_data_filtered, outcome_var, pseudocount = pseudocount)
    test_data_filtered  <- .apply_clr(test_data_filtered, outcome_var, pseudocount = pseudocount)
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    # CLASS BALANCE: computes native per-class weights IF
    # class_balance_method == "class_weights". Kept after the prevalence
    # filter and the CLR, since it only depends on the outcome_var column
    # in train_data_filtered (already filtered), not on the
    # CLR-transformed values of the taxa columns. Unlike the previous
    # approach (case weights via an extra column +
    # workflows::add_case_weights), here no extra column is attached to
    # train_data_filtered -- the weights are passed directly to the
    # engines in .buildWorkflows().
    # ---------------------------------------------------------------
    class_weights_rf <- NULL
    scale_pos_weight_xgb <- NULL
    if (class_balance_method == "class_weights") {
      y <- train_data_filtered[[outcome_var]]
      tab <- table(y)

      # RF (ranger::class.weights): weight inversely proportional to class
      # frequency, named by the factor levels (equivalent to
      # class_weight = "balanced").
      w <- length(y) / (length(tab) * tab)
      class_weights_rf <- stats::setNames(as.numeric(w), names(tab))

      # XGB (xgboost::scale_pos_weight): negative/positive ratio, per
      # xgboost's own default definition. The positive class is always
      # the second level of the outcome factor (see
      # class_of_interest/RUMBLE()).
      negative_class <- levels(y)[1]
      positive_class <- levels(y)[2]
      scale_pos_weight_xgb <- as.numeric(tab[[negative_class]] / tab[[positive_class]])
    }

    # Build recipe and workflows
    set.seed(seed)
    recipe <- .buildRecipe(
      train_data_filtered,
      outcome_var,
      class_balance_method = class_balance_method,
      seed = seed,
      correlation_threshold = correlation_threshold
    )
    workflows <- .buildWorkflows(
      recipe, xgb_trees = xgb_trees, rf_trees = rf_trees, seed = seed,
      class_balance_method = class_balance_method,
      class_weights_rf = class_weights_rf,
      scale_pos_weight_xgb = scale_pos_weight_xgb
    )

    # Tune models (Inner loop)
    tuned <- .tuneModels(
      workflows, train_data_filtered, outcome_var,
      grid_size = grid_size, cv_folds = cv_folds,
      n_cores = n_cores, seed = seed
    )

    # Evaluate on the outer fold test set
    final_models <- .evaluateModels(
      tuned, train_data_filtered, test_data_filtered, outcome_var
    )

    list(
      fold_id = fold_idx,
      train_data = train_data_filtered,
      test_data = test_data_filtered,
      tuned = tuned,
      models = final_models
    )
  })

  ## ==================================================================
  ## 4.5 Aggregate Metrics and Apply Cutoffs
  ## ==================================================================
  msg("\nAggregating metrics across all ", outer_folds, " outer folds...")

  all_metrics <- purrr::map_dfr(fold_results, function(res) {
    purrr::imap_dfr(res$models, function(mod, model_name) {
      mod$metrics %>%
        dplyr::mutate(model = model_name, fold = res$fold_id)
    })
  })

  # Calculate mean and SD of metrics across outer test folds
  metrics_summary <- all_metrics %>%
    dplyr::group_by(model, .metric) %>%
    dplyr::summarize(
      mean_estimate = mean(.estimate, na.rm = TRUE),
      sd_estimate = sd(.estimate, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(Metric = .metric, Mean = mean_estimate, SD = sd_estimate)

  # Extract hyperparameters from each fold for full transparency
  best_hyperparameters <- purrr::map_dfr(fold_results, function(res) {
    purrr::map_dfr(names(res$tuned), function(name) {
      obj <- res$tuned[[name]]
      if (is.null(obj) || is.null(obj$best_params)) return(NULL)
      obj$best_params %>%
        dplyr::select(-dplyr::any_of(".config")) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
        tidyr::pivot_longer(cols = dplyr::everything(),
                            names_to = "Parameter", values_to = "Value") %>%
        dplyr::mutate(Model = name, Fold = res$fold_id) %>%
        dplyr::select(Model, Fold, Parameter, Value)
    })
  })

  # Apply metric_cutoffs based on AVERAGE performance across K-Folds
  models_to_keep <- unique(metrics_summary$model)

  if (!is.null(metric_cutoffs)) {
    msg("Applying metric cutoffs (based on outer folds average)...")
    msg("Cutoff thresholds: ", paste(names(metric_cutoffs), "=", metric_cutoffs, collapse = ", "))

    metrics_for_filter <- metrics_summary %>%
      dplyr::select(model, Metric, Mean) %>%
      tidyr::pivot_wider(names_from = Metric, values_from = Mean)

    models_to_keep <- character(0)
    for (model_name in unique(metrics_summary$model)) {
      model_metrics <- metrics_for_filter[metrics_for_filter$model == model_name, ]
      passes_all_cutoffs <- TRUE

      for (metric_name in names(metric_cutoffs)) {
        if (metric_name %in% colnames(model_metrics)) {
          metric_value <- model_metrics[[metric_name]]
          cutoff_value <- metric_cutoffs[[metric_name]]

          if (is.na(metric_value) || metric_value < cutoff_value) {
            msg("  Model '", model_name, "' failed cutoff for ", metric_name,
                ": ", round(metric_value, 3), " < ", cutoff_value)
            passes_all_cutoffs <- FALSE
            break
          }
        }
      }

      if (passes_all_cutoffs) {
        models_to_keep <- c(models_to_keep, model_name)
        msg("  Model '", model_name, "' passed all cutoffs")
      }
    }

    if (length(models_to_keep) == 0) {
      stop("All models were filtered out by metric_cutoffs. Consider relaxing the thresholds.")
    }
  }

  selected_shap_models <- models_to_keep

  # Validate custom shap_model names against actual trained models
  if (shap_model_is_custom) {
    invalid_models <- setdiff(shap_model, selected_shap_models)
    if (length(invalid_models) > 0) {
      stop("Invalid 'shap_model' selection: ",
           paste(invalid_models, collapse = ", "),
           ". Available models: ", paste(selected_shap_models, collapse = ", "))
    }
  }

  ## ==================================================================
  ## 5. Interpretability (SHAP)
  ## ==================================================================
  # The pipeline always computes SHAP on both the training and test sets:
  compute_train <- TRUE
  compute_test  <- TRUE

  # Consolidation below therefore always builds the union of both sets --
  # the pipeline compares train vs. test SHAP signatures to build the
  # overfitting report.

  msg("[5/5] Computing feature importance (SHAP)...")

  fold_importances <- purrr::map(fold_results, function(res) {
    fold_name <- paste0("Fold_", res$fold_id)
    msg("  -> Processing ", fold_name)

    explainable_models <- res$models[names(res$models) %in% selected_shap_models]
    if (length(explainable_models) == 0) return(NULL)

    # Explainers created on Training data
    explainers <- createExplainers(
      explainable_models, res$train_data, outcome_var,
      class_of_interest = class_of_interest,
      verbose = verbose
    )

    #%% [ORIGINAL -- pre-commit1] Per-fold model_weights removed: it used
    #%% the current fold's bal_accuracy/roc_auc (not the MCC aggregated
    #%% across outer folds), without the w_min documented in Equation 5,
    #%% and it only fed the dead global_importance field of
    #%% computeFeatureImportance() -- it never reached the consensus
    #%% tables reported in the manuscript. The only performance-weighting
    #%% mechanism now is global_weights (Section 5.5 below), which
    #%% faithfully implements Equation 5.
    #%%
    #%% model_weights <- purrr::map_dbl(explainable_models, function(mod) {
    #%%   metrics_df <- mod$metrics
    #%%   bal_acc_row <- metrics_df[metrics_df$.metric == "bal_accuracy", ]
    #%%   if (nrow(bal_acc_row) > 0 && !is.na(bal_acc_row$.estimate[1])) return(bal_acc_row$.estimate[1])
    #%%   auc_row <- metrics_df[metrics_df$.metric == "roc_auc", ]
    #%%   if (nrow(auc_row) > 0 && !is.na(auc_row$.estimate[1])) return(auc_row$.estimate[1])
    #%%   return(1.0)
    #%% })
    #%% model_weights <- pmax(model_weights, 0.01)
    #%% model_weights <- model_weights / sum(model_weights)

    imp_train <- NULL
    imp_test <- NULL

    # SHAP on Training set
    msg("    - Computing SHAP on Train set...")
    imp_train <- computeFeatureImportance(
      explainers, res$train_data, outcome_var,
      class_of_interest = class_of_interest, top_n = top_n,
      repetitions = shap_reps, n_cores = n_cores, shap_method = shap_method,
      verbose = FALSE
    )
    if (!is.null(imp_train) && nrow(imp_train$shap_raw) > 0) {
      imp_train$shap_raw$fold <- res$fold_id
    }

    msg("    - Computing SHAP on Test set...")
    imp_test <- computeFeatureImportance(
      explainers, res$test_data, outcome_var,
      class_of_interest = class_of_interest, top_n = top_n,
      repetitions = shap_reps, n_cores = n_cores, shap_method = shap_method,
      verbose = FALSE
    )
    if (!is.null(imp_test) && nrow(imp_test$shap_raw) > 0) {
      imp_test$shap_raw$fold <- res$fold_id
      imp_test$permutation_raw$fold <- res$fold_id
      imp_test$permutation_top$fold <- res$fold_id
    }

    list(train = imp_train, test = imp_test)
  })

  fold_importances <- fold_importances[!vapply(fold_importances, is.null, logical(1L))]

  ## ==================================================================
  ## 5.5 Consolidate and Summarize Explainability
  ## ==================================================================

  # Consolidate raw SHAP data across folds
  all_shap_raw_train <- purrr::map_dfr(fold_importances, ~ .x$train$shap_raw)
  all_shap_raw_test  <- purrr::map_dfr(fold_importances, ~ .x$test$shap_raw)
  # Determine primary SHAP source for consensus
  primary_shap <- if (compute_test && nrow(all_shap_raw_test) > 0) {
    all_shap_raw_test
  } else {
    all_shap_raw_train
  }

  # SAFETY GATE
  if (nrow(primary_shap) == 0 || !"model" %in% colnames(primary_shap)) {
    msg("  [WARNING] SHAP signatures are empty. Models may have learned nothing.")
    shap_overfitting <- data.frame(model = character(0), spearman_train_test = numeric(0), kendall_train_test = numeric(0))
    model_consistency <- NULL
    selected_shap_models <- character(0)
    feature_frequency <- data.frame()

    importance <- list(
      permutation_raw = data.frame(),
      permutation_top = data.frame(),
      shap_raw_test = data.frame(),
      shap_raw_train = data.frame(),
      shap_top = list(spearman = data.frame()),
      global_importance_test = list(spearman = data.frame()),
      global_importance_train = list(spearman = data.frame())
    )
  } else {

    # Overfitting evaluation (only if both train and test are available)
    if (compute_train && compute_test && nrow(all_shap_raw_train) > 0 && nrow(all_shap_raw_test) > 0) {
      msg("Evaluating SHAP Overfitting (Train vs Test signatures)...")
      shap_overfitting <- evaluateExplainabilityOverfitting(all_shap_raw_train, all_shap_raw_test)
    } else {
      shap_overfitting <- data.frame(model = character(0), spearman_train_test = numeric(0), kendall_train_test = numeric(0))
    }

    # =========================================================================
    # A. MODEL-SPECIFIC SUMMARY (Based on primary SHAP source)
    # =========================================================================
    shap_model_summary <- primary_shap %>%
      dplyr::group_by(model, variable) %>%
      dplyr::summarize(
        mean_abs_contribution = mean(abs(contribution), na.rm = TRUE),
        sd_abs_contribution = sd(abs(contribution), na.rm = TRUE),
        SHAP_Rho = .compute_direction_internal(feature_value, contribution),
        .groups = "drop"
      ) %>%
      dplyr::group_by(model) %>%
      dplyr::arrange(model, dplyr::desc(mean_abs_contribution))

    top_shap_spearman_aggregated <- shap_model_summary %>%
      dplyr::group_by(model) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::ungroup()

    # =========================================================================
    # B. WEIGHTED CONSENSUS SUMMARY
    # =========================================================================
    global_weights <- metrics_summary %>%
      dplyr::filter(Metric == "mcc") %>%
      dplyr::select(model, Mean) %>%
      tibble::deframe()

    if (length(global_weights) == 0 || any(is.na(global_weights))) {
      global_weights <- metrics_summary %>%
        dplyr::filter(Metric == "roc_auc") %>%
        dplyr::select(model, Mean) %>%
        tibble::deframe()
    }
    if (length(global_weights) == 0) {
      global_weights <- stats::setNames(rep(1, length(selected_shap_models)), selected_shap_models)
    }

    global_weights <- pmax(global_weights, 0.01)
    global_weights <- global_weights / sum(global_weights)

    # -------------------------------------------------------------------------
    # B.1 FEATURE SELECTION FREQUENCY (Inter-Folds)
    # -------------------------------------------------------------------------
    msg("Calculating feature selection frequency across folds...")

    feature_frequency <- calculateFeatureFrequency(
      primary_shap, top_n = top_n, n_folds = outer_folds, min_fold_frequency = min_fold_frequency
    )

    # Filter features by min_fold_frequency
    stable_features <- feature_frequency$variable[
      feature_frequency$fold_frequency >= min_fold_frequency
    ]

    msg("  Features passing min_fold_frequency (", min_fold_frequency * 100, "%): ",
        length(stable_features), " / ", nrow(feature_frequency))

    # ---------------------------------------------------------------
    # min_fold_frequency is defined over model x fold evaluations (see
    # calculateFeatureFrequency()), not over outer_folds alone: with M
    # models admitted by metric_cutoffs and outer_folds folds, the only
    # frequencies a feature can actually take are multiples of
    # 1 / (M * outer_folds). A requested min_fold_frequency that is not
    # one of those multiples is silently rounded down in effect (a feature
    # needs ceiling(min_fold_frequency * M * outer_folds) evaluations to
    # pass). This generalizes the fix requested for the specific case of
    # outer_folds = 3, where 0.7 was mistakenly read as "70% of folds"
    # instead of "70% of model x fold evaluations".
    #
    # Gated behind apply_stability_filter: min_fold_frequency only
    # affects the reported results when the filter is actually enabled
    # (see the analogous isTRUE(apply_stability_filter) gate further
    # below, where the filter is applied to primary_global_importance).
    # With apply_stability_filter = FALSE (the default), min_fold_frequency
    # is not used to prune anything, so warning about its reachability
    # unconditionally would print noise for the common case of users who
    # never touch the stability filter at all.
    # ---------------------------------------------------------------
    if (isTRUE(apply_stability_filter) &&
        nrow(feature_frequency) > 0 && "total_possible_evaluations" %in% colnames(feature_frequency)) {
      total_evals <- feature_frequency$total_possible_evaluations[1]
      if (total_evals > 0) {
        step <- 1 / total_evals
        nearest_reachable <- floor((min_fold_frequency + 1e-9) / step) * step
        if (abs(nearest_reachable - min_fold_frequency) > 1e-9) {
          msg("  NOTE: min_fold_frequency = ", min_fold_frequency,
              " is not an exact multiple of 1/", total_evals,
              " (", total_evals, " model x fold evaluations in this run). ",
              "The nearest reachable threshold at or below the requested value is ",
              round(nearest_reachable, 4), " (i.e., a feature must appear in at least ",
              ceiling(min_fold_frequency * total_evals), " of ", total_evals,
              " model x fold evaluations to pass).")
        }
      }
    }

    # -------------------------------------------------------------------------
    # B.2 CONSENSUS FROM TEST (Gold Standard)
    # -------------------------------------------------------------------------
    if (compute_test && nrow(all_shap_raw_test) > 0) {

      # Step 1: Intra-Model Normalization
      all_shap_weighted_test <- all_shap_raw_test %>%
        dplyr::group_by(model) %>%
        dplyr::mutate(
          abs_contribution = abs(contribution),
          # Min-Max Scaling (rescales the absolute SHAP values of all models onto a common 0-1 scale)
          scaled_abs_contribution = (abs_contribution - min(abs_contribution, na.rm = TRUE)) /
            (max(abs_contribution, na.rm = TRUE) - min(abs_contribution, na.rm = TRUE) + 1e-9),
          weight = global_weights[model]
        ) %>%
        dplyr::ungroup()

      # Step 2: Weighted Inter-Model Consensus
      global_importance_spearman_aggregated_test <- all_shap_weighted_test %>%
        dplyr::group_by(variable) %>%
        dplyr::summarise(
          # Uses the scaled variable here!
          mean_abs_contribution = stats::weighted.mean(scaled_abs_contribution, w = weight, na.rm = TRUE),
          #%% [ORIGINAL -- pre-commit1] SHAP_Rho = .compute_direction_internal(feature_value, contribution),
          #%% Unweighted: the direction metric ignored the MCC weight that
          #%% already weighted the magnitude alongside it (line above),
          #%% breaking Equation 5 of the manuscript, which defines a
          #%% single MCC-based weighting scheme applied to the consensus.
          # Computes direction on-the-fly from the original vectors, now
          # weighted by the same `weight` (MCC-based, Eq. 5) used for the magnitude
          SHAP_Rho = .compute_direction_internal(feature_value, contribution, weight = weight),
          n_models  = dplyr::n_distinct(model),
          n_folds   = dplyr::n_distinct(fold),
          .groups   = "drop"
        ) %>%
        dplyr::left_join(
          feature_frequency[, c("variable", "fold_frequency", "n_folds_present", "is_stable")],
          by = "variable"
        ) %>%
        dplyr::arrange(dplyr::desc(mean_abs_contribution))

    } else {
      global_importance_spearman_aggregated_test <- data.frame()
    }

    # -------------------------------------------------------------------------
    # B.3 CONSENSUS FROM TRAIN (Audit / Debugging)
    # -------------------------------------------------------------------------
    if (compute_train && nrow(all_shap_raw_train) > 0) {

      # Step 1: Intra-Model Normalization, replicated for the Training set
      all_shap_weighted_train <- all_shap_raw_train %>%
        dplyr::group_by(model) %>%
        dplyr::mutate(
          abs_contribution = abs(contribution),
          scaled_abs_contribution = (abs_contribution - min(abs_contribution, na.rm = TRUE)) /
            (max(abs_contribution, na.rm = TRUE) - min(abs_contribution, na.rm = TRUE) + 1e-9),
          weight = global_weights[model]
        ) %>%
        dplyr::ungroup()

      # Step 2: Weighted Inter-Model Consensus for the Training set
      global_importance_spearman_aggregated_train <- all_shap_weighted_train %>%
        dplyr::group_by(variable) %>%
        dplyr::summarise(
          # Uses the scaled variable here as well!
          mean_abs_contribution = stats::weighted.mean(scaled_abs_contribution, w = weight, na.rm = TRUE),
          #%% [ORIGINAL -- pre-commit1] SHAP_Rho = .compute_direction_internal(feature_value, contribution),
          #%% Same fix as the test block above: direction is now weighted
          #%% by the MCC-based `weight` (Eq. 5), consistent with the
          #%% magnitude on the previous line.
          SHAP_Rho = .compute_direction_internal(feature_value, contribution, weight = weight),
          n_models  = dplyr::n_distinct(model),
          n_folds   = dplyr::n_distinct(fold),
          .groups   = "drop"
        ) %>%
        dplyr::left_join(
          feature_frequency[, c("variable", "fold_frequency", "n_folds_present", "is_stable")],
          by = "variable"
        ) %>%
        dplyr::arrange(dplyr::desc(mean_abs_contribution))

    } else {
      global_importance_spearman_aggregated_train <- data.frame()
    }

    # =========================================================================
    # FINAL PACKAGING
    # =========================================================================
    importance <- list(
      permutation_raw = if (compute_test) purrr::map_dfr(fold_importances, ~ .x$test$permutation_raw) else data.frame(),

      permutation_top = if (compute_test) {
        purrr::map_dfr(fold_importances, ~ .x$test$permutation_top) %>%
          dplyr::group_by(model, variable) %>%
          # 1. Computes the mean raw delta loss across folds
          dplyr::summarize(mean_delta_loss = mean(mean_delta_loss, na.rm = TRUE), .groups = "drop") %>%
          # 2. Groups by model to apply Min-Max Scaling
          dplyr::group_by(model) %>%
          dplyr::mutate(
            # Normalizes the scale (0 to 1) to ensure a fair inter-model comparison
            mean_delta_loss = (mean_delta_loss - min(mean_delta_loss, na.rm = TRUE)) /
              (max(mean_delta_loss, na.rm = TRUE) - min(mean_delta_loss, na.rm = TRUE) + 1e-9)
          ) %>%
          dplyr::ungroup() %>%
          dplyr::arrange(model, dplyr::desc(mean_delta_loss))
      } else { data.frame() },

      shap_raw_test = all_shap_raw_test,
      shap_raw_train = all_shap_raw_train,
      shap_top = list(spearman = top_shap_spearman_aggregated),
      global_importance_test = list(spearman = global_importance_spearman_aggregated_test),
      global_importance_train = list(spearman = global_importance_spearman_aggregated_train)
    )
  }

  # ==================================================================
  # SHAP CONSISTENCY (Inter-Folds)
  # ==================================================================
  msg("Calculating SHAP signature consistency across folds...")

  shap_consistency <- if (nrow(primary_shap) > 0 && "fold" %in% colnames(primary_shap)) {
    .calculateShapConsistency(primary_shap, top_n = top_n)
  } else {
    data.frame(model = character(0), mean_spearman = numeric(0), mean_jaccard_top20 = numeric(0))
  }

  # C. CONSOLIDATE MODELS FOR PLOTTING (Out-of-Fold)
  final_models_consolidated <- purrr::map(
    stats::setNames(unique(all_metrics$model), unique(all_metrics$model)),
    function(model_name) {
      combined_preds <- purrr::map_dfr(fold_results, function(res) {
        if (model_name %in% names(res$models)) {
          preds <- res$models[[model_name]]$predictions
          preds$fold <- res$fold_id
          return(preds)
        }
        return(NULL)
      })
      formatted_metrics <- metrics_summary %>%
        dplyr::filter(model == model_name) %>%
        dplyr::transmute(.metric = Metric, .estimate = Mean)

      list(predictions = combined_preds, metrics = formatted_metrics)
    }
  )

  # Inner loop CV metrics
  cv_metrics_summary <- purrr::map_dfr(fold_results, function(res) {
    purrr::map_dfr(names(res$tuned), function(name) {
      obj <- res$tuned[[name]]
      if (is.null(obj) || !inherits(obj$tuning_results, "tune_results")) return(NULL)

      # 1. Isolates the ID of the configuration that actually won the tuning
      best_config_id <- obj$best_params$.config[1]

      tune::collect_metrics(obj$tuning_results) %>%
        # 2. STRICT FILTER: keeps only the metrics tied to the winning hyperparameters
        dplyr::filter(.config == best_config_id) %>%
        dplyr::select(.metric, mean, std_err) %>%
        dplyr::rename(Metric = .metric, CV_Mean = mean, CV_StdErr = std_err) %>%
        dplyr::mutate(Model = name, Fold = res$fold_id) %>%
        dplyr::select(Model, Fold, Metric, CV_Mean, CV_StdErr)
    })
  })

  # Inter-Model Consistency
  msg("Evaluating Inter-Model Consistency (Jaccard & Spearman)...")
  model_consistency <- if (nrow(primary_shap) > 0) {
    evaluateInterModelConsistency(primary_shap, top_n_jaccard = top_n)
  } else {
    NULL
  }

  # Hyperparameter Stability Summary
  hyperparameter_stability <- .summarizeHyperparameterStability(best_hyperparameters)

  ## ==================================================================
  ## Plot and Table Generation
  ## ==================================================================
  msg("Generating plots based on consolidated out-of-fold data...")

  plot_filtered_counts <- filtered_counts
  plot_metadata <- prep$metadata

  # Determine primary importance for plotting
  primary_global_importance <- if (nrow(importance$global_importance_test$spearman) > 0) {
    importance$global_importance_test$spearman
  } else {
    importance$global_importance_train$spearman
  }

  if (isTRUE(apply_stability_filter) && "is_stable" %in% colnames(primary_global_importance)) {
    n_before <- nrow(primary_global_importance)
    primary_global_importance <- primary_global_importance %>%
      dplyr::filter(is_stable)
    msg("Stability filter ON: retained ", nrow(primary_global_importance),
        " / ", n_before, " features (fold_frequency >= ", min_fold_frequency, ")")
  }

  # Integrated biomarker table (primary)
  integrated_table_consensus <- if (nrow(primary_global_importance) > 0) {
    getIntegratedBiomarkerTable(
      filtered_counts = plot_filtered_counts,
      metadata = plot_metadata,
      global_importance = primary_global_importance,
      target_var = outcome_var,
      top_n = top_n,
      class_of_interest = class_of_interest,
      negative_class = negative_class,
      da_results = da_results,
      metric_name = "Spearman",
      min_abundance = min_abundance
    )
  } else {
    data.frame()
  }

  # Model-specific plots
  if (length(shap_model) == 1L && identical(shap_model, "all")) {
    shap_models_to_plot <- selected_shap_models
  } else if (length(shap_model) == 1L && identical(shap_model, "consensus")) {
    shap_models_to_plot <- character(0)
  } else {
    shap_models_to_plot <- intersect(shap_model, selected_shap_models)
  }

  model_specific_plots <- list(
    shap_spearman = setNames(vector("list", length(shap_models_to_plot)), shap_models_to_plot),
    shap_beeswarm = setNames(vector("list", length(shap_models_to_plot)), shap_models_to_plot),
    shap_dependence = setNames(vector("list", length(shap_models_to_plot)), shap_models_to_plot),
    biomarker_integrated_spearman = setNames(vector("list", length(shap_models_to_plot)), shap_models_to_plot)
  )

  if (length(shap_models_to_plot) > 0L && nrow(primary_global_importance) > 0) {
    for (model_name in shap_models_to_plot) {
      model_specific_plots$shap_spearman[[model_name]] <- plotShapGlobal(
        importance$shap_top$spearman, outcome_var,
        class_of_interest = class_of_interest, negative_class = negative_class,
        top_n = top_n, metric_name = "Spearman", model = model_name
      )
      model_specific_plots$shap_beeswarm[[model_name]] <- plotShapBeeswarm(
        primary_shap, prediction_data = NULL, target_var = outcome_var,
        top_n = top_n, model = model_name
      )
      model_specific_plots$shap_dependence[[model_name]] <- plotShapDependence(
        primary_shap, top_n = top_n, model = model_name
      )
      model_specific_plots$biomarker_integrated_spearman[[model_name]] <- plotBiomarkerIntegrated(
        filtered_counts = plot_filtered_counts, metadata = plot_metadata,
        global_importance = importance$shap_top$spearman, target_var = outcome_var,
        top_n = top_n, class_of_interest = class_of_interest, negative_class = negative_class,
        da_results = da_results, metric_name = "Spearman", min_abundance = min_abundance, model = model_name
      )
    }
  }

  # Build main plots list
  plots <- list(
    metrics   = plotMetricsComparison(final_models_consolidated, fold_metrics = all_metrics),
    roc       = plotRocCurves(final_models_consolidated, outcome_var),
    cm        = plotConfusionMatrices(final_models_consolidated, outcome_var),
    oof_density = plotOOFDensity(final_models_consolidated, outcome_var),
    shap_spearman = if (nrow(primary_global_importance) > 0) {
      plotShapGlobal(
        primary_global_importance, outcome_var,
        class_of_interest = class_of_interest, negative_class = negative_class,
        top_n = top_n, metric_name = "Spearman"
      )
    } else { NULL },
    shap_beeswarm = if (nrow(primary_shap) > 0) {
      plotShapBeeswarm(
        primary_shap, prediction_data = NULL, target_var = outcome_var, top_n = top_n
      )
    } else { NULL },
    shap_dependence = if (nrow(primary_shap) > 0) {
      plotShapDependence(primary_shap, top_n = top_n)
    } else { NULL },
    shap_prevalence = if (nrow(primary_global_importance) > 0) {
      plotTaxaPrevalence(
        filtered_counts = plot_filtered_counts, metadata = plot_metadata,
        global_importance = primary_global_importance, target_var = outcome_var,
        min_abundance = min_abundance, top_n = top_n
      )
    } else { NULL },
    biomarker_integrated_spearman = if (nrow(primary_global_importance) > 0) {
      plotBiomarkerIntegrated(
        filtered_counts = plot_filtered_counts, metadata = plot_metadata,
        global_importance = primary_global_importance, target_var = outcome_var,
        top_n = top_n, class_of_interest = class_of_interest, negative_class = negative_class,
        da_results = da_results, metric_name = "Spearman", min_abundance = min_abundance
      )
    } else { NULL },
    shap_vs_permutation = if (nrow(primary_global_importance) > 0 && nrow(importance$permutation_top) > 0) {
      plotSHAPvsPermutation(primary_global_importance, importance$permutation_top, feature_frequency)
    } else { NULL },
    model_specific = model_specific_plots
  )

  metrics_summary_table <- summariseMetrics(final_models_consolidated)

  ## ==================================================================
  ## Export Files
  ## ==================================================================
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    prefix <- file.path(output_dir, project_name)
    msg("Saving outputs to: ", output_dir)

    plot_names <- c("metrics", "roc", "cm", "oof_density",
                    "shap_spearman", "shap_beeswarm", "shap_dependence",
                    "shap_prevalence", "biomarker_integrated_spearman",
                    "shap_vs_permutation")

    plot_sizes <- list(
      metrics = c(8, 5), roc = c(8, 5), cm = c(10, 8),
      oof_density = c(10, 6), shap_spearman = c(9, 8),
      shap_beeswarm = c(10, 8), shap_dependence = c(12, 10),
      shap_prevalence = c(9, 8),
      biomarker_integrated_spearman = c(12, 9),
      shap_vs_permutation = c(9, 7)
    )

    for (pname in plot_names) {
      if (!is.null(plots[[pname]])) {
        sz <- plot_sizes[[pname]]
        ggplot2::ggsave(
          filename = paste0(prefix, "_", pname, ".png"),
          plot = plots[[pname]], width = sz[1L], height = sz[2L], dpi = 300
        )
      }
    }

    if (length(shap_models_to_plot) > 0L) {
      for (model_name in shap_models_to_plot) {
        model_slug <- gsub("[^A-Za-z0-9]+", "_", model_name)

        # 1. Saves the model's global SHAP bar plot
        if (!is.null(plots$model_specific$shap_spearman[[model_name]])) {
          ggplot2::ggsave(filename = paste0(prefix, "_", model_slug, "_shap_spearman.png"), plot = plots$model_specific$shap_spearman[[model_name]], width = 9, height = 8, dpi = 300)
        }

        # 2. Saves the model's beeswarm plot
        if (!is.null(plots$model_specific$shap_beeswarm[[model_name]])) {
          ggplot2::ggsave(filename = paste0(prefix, "_", model_slug, "_shap_beeswarm.png"), plot = plots$model_specific$shap_beeswarm[[model_name]], width = 10, height = 8, dpi = 300)
        }

        # 3. Saves the model's SHAP dependence plot
        if (!is.null(plots$model_specific$shap_dependence[[model_name]])) {
          ggplot2::ggsave(filename = paste0(prefix, "_", model_slug, "_shap_dependence.png"), plot = plots$model_specific$shap_dependence[[model_name]], width = 12, height = 10, dpi = 300)
        }

        # 4. FIX: saves the model's Integrated Biomarker Dashboard
        if (!is.null(plots$model_specific$biomarker_integrated_spearman[[model_name]])) {
          ggplot2::ggsave(filename = paste0(prefix, "_", model_slug, "_biomarker_integrated_spearman.png"), plot = plots$model_specific$biomarker_integrated_spearman[[model_name]], width = 12, height = 9, dpi = 300)
        }
      }
    }

    # =========================================================================
    # TSV EXPORTS
    # =========================================================================

    # 1. Sample-level Out-of-Fold (OOF) predictions
    msg("Exporting Out-of-Fold predictions...")
    oof_preds <- purrr::imap_dfr(final_models_consolidated, function(mod_data, model_name) {
      df <- mod_data$predictions
      df$model <- model_name
      return(df)
    })

    # Reordering columns for a cleaner, more readable layout
    oof_preds <- oof_preds %>%
      dplyr::select(sample_id, fold, model, !!rlang::sym(outcome_var), .pred_class, dplyr::starts_with(".pred_"))

    readr::write_tsv(oof_preds, paste0(prefix, "_oof_predictions.tsv"))
    # --- NEW BLOCK: AUTOMATED PREDICTIVE DIAGNOSTICS ---
    oof_diagnostics <- analyzeOOFPredictions(oof_preds, outcome_var)
    if (!is.null(oof_diagnostics)) {
      readr::write_tsv(oof_diagnostics$confidence, paste0(prefix, "_oof_confidence_summary.tsv"))
      if (nrow(oof_diagnostics$hard_samples) > 0) {
        readr::write_tsv(oof_diagnostics$hard_samples, paste0(prefix, "_oof_hard_samples.tsv"))
      }
    }
    # -------------------------------------------------------
    readr::write_tsv(metrics_summary_table, paste0(prefix, "_model_metrics.tsv"))
    readr::write_tsv(best_hyperparameters, paste0(prefix, "_best_hyperparameters.tsv"))
    readr::write_tsv(cv_metrics_summary, paste0(prefix, "_cv_metrics.tsv"))
    readr::write_tsv(all_metrics, paste0(prefix, "_model_metrics_per_fold.tsv"))

    # Unified Tabular Export: SHAP vs. Permutation
    if (nrow(primary_global_importance) > 0 && nrow(importance$permutation_top) > 0) {
      shap_vs_perm_df <- dplyr::inner_join(
        primary_global_importance %>% dplyr::rename(shap_importance = mean_abs_contribution),
        importance$permutation_top %>%
          dplyr::group_by(variable) %>%
          dplyr::summarize(perm_importance = mean(mean_delta_loss, na.rm = TRUE), .groups = "drop"),
        by = "variable"
      ) %>%
        dplyr::arrange(dplyr::desc(shap_importance), dplyr::desc(perm_importance))

      readr::write_tsv(shap_vs_perm_df, paste0(prefix, "_shap_vs_permutation.tsv"))
    }

    # 2. Hyperparameter Stability
    if (nrow(hyperparameter_stability) > 0) {
      readr::write_tsv(hyperparameter_stability, paste0(prefix, "_hyperparameter_stability.tsv"))
    }

    # 3. SHAP & Biomarkers (primary source)
    if (nrow(primary_global_importance) > 0) {
      readr::write_tsv(primary_global_importance, paste0(prefix, "_shap_global_consensus.tsv"))
      readr::write_tsv(integrated_table_consensus, paste0(prefix, "_integrated_biomarkers.tsv"))
    }
    if (nrow(importance$shap_top$spearman) > 0) {
      readr::write_tsv(importance$shap_top$spearman, paste0(prefix, "_shap_model_specific.tsv"))
    }
    if (nrow(importance$permutation_top) > 0) {
      readr::write_tsv(importance$permutation_top, paste0(prefix, "_permutation_top.tsv"))
    }

    # 4. Raw SHAP data
    if (nrow(all_shap_raw_test) > 0) {
      readr::write_tsv(all_shap_raw_test, paste0(prefix, "_shap_raw_TEST.tsv"))
    }
    if (nrow(all_shap_raw_train) > 0) {
      readr::write_tsv(all_shap_raw_train, paste0(prefix, "_shap_raw_TRAIN.tsv"))
    }
    # --- NEW BLOCK: SHAP MAGNITUDE STABILITY (STEP 4) ---
    if (nrow(all_shap_raw_train) > 0 && nrow(all_shap_raw_test) > 0) {
      shap_magnitude_stab <- evaluateShapMagnitudeStability(
        shap_raw_train = all_shap_raw_train,
        shap_raw_test = all_shap_raw_test,
        top_n = top_n
      )
      if (!is.null(shap_magnitude_stab)) {
        readr::write_tsv(shap_magnitude_stab, paste0(prefix, "_shap_magnitude_stability_top", top_n, ".tsv"))
      }
    }
    # ------------------------------------------------------------

    # 5. Consistency & Stability Metrics
    if (nrow(shap_overfitting) > 0) {
      readr::write_tsv(shap_overfitting, paste0(prefix, "_shap_overfitting_train_vs_test.tsv"))
    }
    if (nrow(shap_consistency) > 0) {
      readr::write_tsv(shap_consistency, paste0(prefix, "_shap_stability_across_folds.tsv"))
    }
    if (!is.null(model_consistency)) {
      readr::write_tsv(model_consistency, paste0(prefix, "_inter_model_consistency.tsv"))
    }
    if (nrow(feature_frequency) > 0) {
      readr::write_tsv(feature_frequency, paste0(prefix, "_feature_selection_frequency.tsv"))
    }

    # 6. Differential Abundance
    if (run_da && !is.null(da_results) && nrow(integrated_table_consensus) > 0) {
      top_da_results <- da_results %>%
        dplyr::filter(Taxon %in% integrated_table_consensus$Taxon)
      readr::write_tsv(top_da_results, paste0(prefix, "_top_differential_abundance.tsv"))

      # Full DA table (mandatory supplementary file for submission)
      readr::write_tsv(da_results, paste0(prefix, "_full_differential_abundance.tsv"))
    }

    # 7. Cleaned and Aligned Data Matrices (Data Availability / Reproducibility)
    msg("Exporting cleaned data matrices (Data Availability)...")

    # 7.1 Filtered abundance matrix (Raw/Relative)
    readr::write_tsv(plot_filtered_counts %>% tibble::rownames_to_column("sample_id"),
                     paste0(prefix, "_filtered_counts.tsv"))

    # 7.2 Final aligned metadata (excluding samples that had NA in the outcome)
    readr::write_tsv(plot_metadata %>% tibble::rownames_to_column("sample_id"),
                     paste0(prefix, "_aligned_metadata.tsv"))

    # 7.3 Global CLR matrix (computed without leakage, for use in PCoA, PERMANOVA, or network analyses in Python)
    global_clr <- .apply_clr(analysis_df, outcome_var, pseudocount = pseudocount) %>%
      dplyr::select(-!!rlang::sym(outcome_var)) %>%
      tibble::rownames_to_column("sample_id")
    readr::write_tsv(global_clr, paste0(prefix, "_global_clr_features.tsv"))

    # Build the publication-formatted tables
    outer_summary <- summarizeOuterCV(all_metrics)
    inner_summary <- summarizeInnerCV(cv_metrics_summary)

    if (!is.null(outer_summary)) {
      readr::write_tsv(outer_summary$publication_table, paste0(prefix, "_publication_outer_cv_consistency.tsv"))
      readr::write_tsv(outer_summary$raw_stats, paste0(prefix, "_raw_outer_cv_discrepancy.tsv"))
    }

    if (!is.null(inner_summary)) {
      readr::write_tsv(inner_summary$publication_table, paste0(prefix, "_publication_inner_cv_consistency.tsv"))
    }
    # --- NEW BLOCK: THE GENERALIZATION GAP ---
    if (!is.null(outer_summary) && !is.null(inner_summary)) {
      gap_summary <- evaluateGeneralizationGap(
        outer_raw = outer_summary$raw_stats,
        inner_raw = inner_summary$raw_stats
      )
      if (!is.null(gap_summary)) {
        msg("Exporting Generalization Gap (Predictive Overfitting) table...")
        readr::write_tsv(gap_summary, paste0(prefix, "_publication_generalization_gap.tsv"))
      }
    }
    # --- NEW BLOCK: SHAP RANK STABILITY (Top N) ---
    if (nrow(all_shap_raw_train) > 0 && nrow(all_shap_raw_test) > 0) {
      shap_rank_stab <- evaluateShapRankStability(
        shap_raw_train = all_shap_raw_train,
        shap_raw_test = all_shap_raw_test,
        top_n = top_n
      )
      if (!is.null(shap_rank_stab)) {
        readr::write_tsv(shap_rank_stab, paste0(prefix, "_shap_rank_stability_top", top_n, ".tsv"))
      }
    }
    # ------------------------------------------------------------
    # -----------------------------------------------------------
    # --- NEW BLOCK: CORE CONSENSUS BIOMARKERS (STEP 5) ---
    if (!is.null(importance$shap_top$spearman) && nrow(importance$shap_top$spearman) > 0) {
      core_consensus_res <- extractCoreConsensus(importance$shap_top$spearman, top_n = top_n)

      if (!is.null(core_consensus_res) && nrow(core_consensus_res$core_consensus) > 0) {
        msg("Exporting Core Consensus Biomarkers (Intersect across ", core_consensus_res$n_models, " models)...")
        readr::write_tsv(core_consensus_res$core_consensus, paste0(prefix, "_core_consensus_biomarkers.tsv"))
      } else {
        msg("No core consensus features found across all models in the top ", top_n, ".")
      }
    }
    # ------------------------------------------------------------

    # 8. Trained Model Objects and R Structures (Serialization)
    msg("Serializing trained models and objects to RDS...")

    # 8.1 Saves the consolidated outputs (OOF predictions, summary metrics, and generated plots)
    pipeline_outputs <- list(
      models = final_models_consolidated,
      metrics = metrics_summary,
      plots = plots,
      integrated_table = integrated_table_consensus
    )
    saveRDS(pipeline_outputs, paste0(prefix, "_pipeline_outputs.rds"))

    # 8.2 Saves the object holding the actual fitted models per fold (essential for inference/external validation)
    saveRDS(fold_results, paste0(prefix, "_trained_fold_models.rds"))

    msg("All outputs saved.")
  }

  msg("======================================================")
  msg("RUMBLE Pipeline complete (Nested CV Mode).")
  msg("======================================================")

  # Final return value of the RUMBLE function
  list(
    models              = final_models_consolidated,
    metrics             = metrics_summary,
    importance          = importance,
    cv_metrics          = cv_metrics_summary,
    hyperparameters     = best_hyperparameters,
    plots               = plots,
    selected_models     = selected_shap_models,
    shap_overfitting    = shap_overfitting,
    model_consistency   = model_consistency,
    da_results          = da_results,
    integrated_table    = integrated_table_consensus,
    shap_consistency    = shap_consistency,
    feature_frequency   = feature_frequency,
    hyperparameter_stability = hyperparameter_stability,
    metrics_per_fold    = all_metrics
  )
}

## ==================================================================
## Internal Helper Functions
## ==================================================================

# REMOVED: .compute_direction (replaced by .compute_direction_internal in explain_functions.R)


#' Calculate SHAP Consistency Across Folds
#'
#' Computes pairwise Spearman correlation and Jaccard similarity of SHAP
#' rankings between all pairs of outer folds for each model.
#'
#' @param shap_raw Data.frame of raw SHAP values with fold column.
#' @param top_n Integer. Number of top features for Jaccard.
#' @return Data.frame with model, mean_spearman, mean_jaccard_top20.
#' @noRd
.calculateShapConsistency <- function(shap_raw, top_n = 20L) {
  purrr::map_dfr(unique(shap_raw$model), function(mod_name) {
    df_mod <- shap_raw[shap_raw$model == mod_name, ]
    folds <- unique(df_mod$fold)

    if (length(folds) < 2) {
      return(data.frame(model = mod_name, mean_spearman = NA_real_, mean_jaccard_top20 = NA_real_))
    }

    fold_signatures <- purrr::map(folds, function(f) {
      df_f <- df_mod[df_mod$fold == f, ]
      ranked <- df_f %>%
        dplyr::group_by(variable) %>%
        dplyr::summarize(importance = mean(abs(contribution), na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(importance))

      list(
        fold = f,
        top_n = utils::head(ranked$variable, top_n),
        scores = stats::setNames(ranked$importance, ranked$variable)
      )
    })

    pairs <- utils::combn(length(folds), 2)

    metrics_pair <- purrr::map_dfr(seq_len(ncol(pairs)), function(i) {
      sig1 <- fold_signatures[[pairs[1, i]]]
      sig2 <- fold_signatures[[pairs[2, i]]]

      # Jaccard of Top-N
      intersect_len <- length(intersect(sig1$top_n, sig2$top_n))
      union_len <- length(union(sig1$top_n, sig2$top_n))
      jaccard <- intersect_len / union_len

      # Spearman of full ranking
      common_vars <- intersect(names(sig1$scores), names(sig2$scores))
      rho <- if (length(common_vars) > 3) {
        suppressWarnings(stats::cor(sig1$scores[common_vars], sig2$scores[common_vars], method = "spearman"))
      } else {
        NA_real_
      }

      data.frame(jaccard = jaccard, spearman = rho)
    })

    data.frame(
      model = mod_name,
      mean_spearman = mean(metrics_pair$spearman, na.rm = TRUE),
      mean_jaccard_top20 = mean(metrics_pair$jaccard, na.rm = TRUE)
    )
  })
}

#' Summarize Hyperparameter Stability Across Folds
#'
#' Calculates the mode, variance, and coefficient of variation for each
#' hyperparameter across outer folds.
#'
#' @param best_hyperparameters Data.frame with columns Model, Fold, Parameter, Value.
#' @return Data.frame with stability metrics.
#' @noRd
.summarizeHyperparameterStability <- function(best_hyperparameters) {
  if (nrow(best_hyperparameters) == 0) return(data.frame())

  best_hyperparameters %>%
    dplyr::group_by(Model, Parameter) %>%
    dplyr::summarize(
      n_folds = dplyr::n(),
      n_unique_values = dplyr::n_distinct(Value),
      mode_value = names(sort(table(Value), decreasing = TRUE))[1],
      values_list = paste(Value, collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      stability_flag = dplyr::case_when(
        n_unique_values == 1 ~ "Stable",
        n_unique_values <= 2 ~ "Moderate",
        TRUE ~ "Unstable"
      )
    ) %>%
    dplyr::arrange(Model, dplyr::desc(n_unique_values))
}
#' Apply Centered Log-Ratio (CLR) Transformation
#'
#' @description Internal function to apply CLR intra-fold to avoid data leakage.
#' @param data Data.frame containing features and target variable.
#' @param target_var Character name of the target variable.
#' @param pseudocount Numeric pseudocount to handle exact zeros.
#'
#' @keywords internal
#' @return A data.frame with CLR-transformed features.
.apply_clr <- function(data, target_var, pseudocount = 1e-6) {
  taxa_cols <- setdiff(colnames(data), target_var)
  mat <- as.matrix(data[, taxa_cols])
  mat[mat == 0] <- pseudocount
  log_mat <- log(mat)
  clr_mat <- sweep(log_mat, 1, rowMeans(log_mat), "-")
  data[, taxa_cols] <- clr_mat
  return(data)
}
