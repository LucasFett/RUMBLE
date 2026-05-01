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
#' @param shap_data Character. Which data to use for SHAP calculation and prevalence plotting.
#'   This parameter is mandatory and must be explicitly specified by the user.
#'   Options are \code{"train"} (recommended for biomarker discovery, uses training data),
#'   \code{"test"} (for prediction-focused analysis), or \code{"all"} (uses entire dataset).
#'   The choice reflects the analytical objective: use \code{"train"} to map what the model
#'   learned during training, \code{"test"} for held-out performance evaluation, or
#'   \code{"all"} for maximum statistical power in exploratory analysis.
#' @param shap_reps Integer. Number of repetitions for SHAP and
#'   permutation importance calculations (default 25).
#'   Higher values provide more stable feature importance rankings.
#'   The default value (25) offers an excellent balance between robustness
#'   and computation time. For critical studies or meta-analyses, consider 100
#'   or more repetitions.
#'   \strong{Note:} This parameter significantly affects computation time.
#' @param shap_method Character. Method for SHAP calculation. Options are
#'   \code{"exact"} (default, uses DALEX for rigorous SHAP values) or
#'   \code{"fast"} (uses fastshap package for faster approximation).
#'   The exact method is recommended for publication-quality results.
#'   The fast method provides a good approximation with significantly
#'   reduced computation time, suitable for exploratory analysis.
#' @param shap_model Character vector controlling which model-specific SHAP
#'   profiles are generated in addition to the consensus plots. Use
#'   \code{"all"} (default) to generate isolated outputs for every model that
#'   passes \code{metric_cutoffs}, \code{"consensus"} to skip isolated model
#'   plots, or provide one or more model names (for example,
#'   \code{c("RF", "XGB")}).
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
#' @param normalization_method Character. Method for compositional normalization.
#'   Currently only \code{"clr"} (Centered Log-Ratio) is supported (default).
#'   The CLR transformation is recommended for all machine learning
#'   algorithms (Random Forest, XGBoost, KNN, ENET) as it preserves mathematical monotonicity
#'   and ensures proper feature importance calculations via SHAP values.
#' @param seed Integer. Random seed for reproducibility (default 42).
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#' @param metric_cutoffs Named numeric vector. Optional quality filter for models
#'   based on performance metrics. Models failing to meet the specified thresholds
#'   are excluded from the SHAP analysis. Example:
#'   \code{metric_cutoffs = c("roc_auc" = 0.75, "f_meas" = 0.60)}.
#'   Default is \code{NULL} (no filtering). This parameter ensures that only
#'   high-quality models contribute to the consensus biomarker ranking.
#' @param run_da Logical. Whether to run differential abundance analysis
#'   (default \code{TRUE}). When \code{TRUE}, the integrated biomarker plots
#'   will include a third panel with differential abundance metrics (Log2FC).
#' @param da_method Character. Differential abundance method to use.
#'   Options are \code{"wilcoxon"} (default fallback requiring no Bioconductor dependencies),
#'   \code{"ancombc"}, or \code{"aldex2"}. Only used if \code{run_da = TRUE}.
#'
#' @return A named list with seven elements:
#'   \itemize{
#'     \item \code{models}: Named list of fitted model objects and predictions.
#'     \item \code{metrics}: Data.frame summarising performance metrics.
#'     \item \code{importance}: List containing SHAP and permutation results.
#'     \item \code{plots}: List of ggplot objects, including consensus and
#'       optional model-specific SHAP plots.
#'     \item \code{data}: List with train and test data.frames.
#'     \item \code{selected_models}: Character vector of models used for
#'       interpretability after applying \code{metric_cutoffs}.
#'     \item \code{da_results}: Data.frame with differential abundance results
#'       (if \code{run_da = TRUE}), or \code{NULL} otherwise.
#'   }
#' @examples
#'   # Load example data provided by the package
#'   ps_path <- system.file("extdata",
#'                          "PRJEB38465_phyloseq_com_metadados_completos.rds",
#'                          package = "RUMBLE")
#'
#'   if (file.exists(ps_path)) {
#'     ps <- readRDS(ps_path)
#'
#'     # EXAMPLE 1: Run pipeline with minimal parameters for demonstration speed
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
#'       verbose = FALSE,
#'       shap_data = "train",
#'       xgb_trees = 5,
#'       rf_trees = 5
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
#'       shap_data = "train",
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
                   shap_reps = 25L,
                   xgb_trees = 1000L,
                   rf_trees = 500L,
                   min_prevalence = 0.05,
                   min_abundance = 0.0001,
                   remove_unclassified = FALSE,
                   normalization_method = "clr",
                   seed = 42L,
                   verbose = TRUE,
                   shap_data,
                   shap_method = "exact",
                   shap_model = "all",
                   metric_cutoffs = NULL,
                   run_da = TRUE,
                   da_method = c("wilcoxon", "ancombc", "aldex2")) {

  msg <- function(...) {
    if (verbose) message(...)
  }

  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  if (missing(shap_data)) {
    stop("Parameter 'shap_data' is mandatory. Please specify one of: ",
         "'train' (recommended for biomarker discovery), ",
         "'test' (for prediction-focused analysis), or ",
         "'all' (for maximum statistical power).")
  }
  # Validate da_method parameter
  da_method <- match.arg(da_method)

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
    normalization_method = normalization_method,
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

  ## Remove samples with NA in outcome (CRITICAL FIX: Synchronize all dataframes)
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
  # This ensures that SHAP values are calculated with respect to the class of interest
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

    # Run DA analysis (directly using the properly formatted prep$metadata)
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
  set.seed(seed)
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
  set.seed(seed)
  final_models <- .evaluateModels(
    tuned, train_data, test_data, outcome_var
  )

  ## ==================================================================
  ## 4.5 Apply metric cutoffs (optional quality filter)
  ## ==================================================================

  # Por padrão, todos os modelos vão para a explicabilidade
  explainable_models <- final_models

  if (!is.null(metric_cutoffs)) {
    msg("Applying metric cutoffs for model selection...")
    msg("Cutoff thresholds: ", paste(names(metric_cutoffs), "=", metric_cutoffs, collapse = ", "))

    # Build a summary of metrics for filtering
    metrics_for_filter <- purrr::imap_dfr(final_models, function(obj, name) {
      obj$metrics %>%
        dplyr::select(.metric, .estimate) %>%
        tidyr::pivot_wider(names_from = .metric, values_from = .estimate) %>%
        dplyr::mutate(model = name)
    })

    # Check which models pass all cutoffs
    models_to_keep <- character(0)
    for (model_name in names(final_models)) {
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
        } else {
          warning("Metric '", metric_name, "' not found in model results. Skipping this cutoff.")
        }
      }

      if (passes_all_cutoffs) {
        models_to_keep <- c(models_to_keep, model_name)
        msg("  Model '", model_name, "' passed all cutoffs")
      }
    }

    # Filter explainable_models (mantendo final_models intacto)
    n_models_before <- length(final_models)
    explainable_models <- final_models[models_to_keep]
    n_models_after <- length(explainable_models)

    if (n_models_after == 0) {
      stop("All models were filtered out by metric_cutoffs. ",
           "Consider relaxing the thresholds.")
    }

    msg("Models retained after filtering (for SHAP): ", n_models_after, " out of ", n_models_before)
  }

  ## ==================================================================
  ## 5. Interpretability
  ## ==================================================================
  msg("[5/5] Computing feature importance...")

  # Determine which data to use for SHAP predictions
  if (shap_data == "all") {
    prediction_data <- analysis_df
    msg("Using ALL samples for SHAP interpretation (default, maximum power)")
  } else if (shap_data == "test") {
    prediction_data <- test_data
    msg("Using TEST set for SHAP interpretation")
  } else if (shap_data == "train") {
    prediction_data <- train_data
    msg("Using TRAIN set for SHAP interpretation")
  } else {
    stop("Invalid 'shap_data' parameter. Use 'all', 'test', or 'train'.")
  }

  selected_shap_models <- names(explainable_models)

  if (length(selected_shap_models) == 0L) {
    stop("No models available for interpretability after filtering.")
  }

  if (length(shap_model) == 1L && identical(shap_model, "all")) {
    shap_models_to_plot <- selected_shap_models
  } else if (length(shap_model) == 1L && identical(shap_model, "consensus")) {
    shap_models_to_plot <- character(0)
  } else {
    invalid_models <- setdiff(shap_model, selected_shap_models)
    if (length(invalid_models) > 0L) {
      stop(
        "Invalid 'shap_model' selection: ",
        paste(invalid_models, collapse = ", "),
        ". Available models after filtering: ",
        paste(selected_shap_models, collapse = ", ")
      )
    }
    shap_models_to_plot <- shap_model
  }

  # Explainer MUST always be created with train data
  # Utilizando a lista explainable_models filtrada
  explainers <- createExplainers(
    explainable_models, train_data, outcome_var,
    class_of_interest = class_of_interest
  )

  importance <- computeFeatureImportance(
    explainers, prediction_data, outcome_var,
    class_of_interest = class_of_interest,
    top_n = top_n, repetitions = shap_reps,
    n_cores = n_cores, shap_method = shap_method,
    verbose = verbose
  )

  ## ==================================================================
  ## Generate plots
  ## ==================================================================
  msg("Generating plots...")

  # Sincroniza as amostras dos gráficos com as amostras usadas no SHAP
  plot_samples <- rownames(prediction_data)
  plot_filtered_counts <- filtered_counts[plot_samples, , drop = FALSE]
  plot_metadata <- prep$metadata[plot_samples, , drop = FALSE]

  msg("Generating integrated data tables...")
  # Extração da Tabela Integrada (Consensus)
  integrated_table_consensus <- getIntegratedBiomarkerTable(
    filtered_counts = plot_filtered_counts,
    metadata = plot_metadata,
    global_importance = importance$global_importance$spearman,
    target_var = outcome_var,
    top_n = top_n,
    class_of_interest = class_of_interest,
    negative_class = negative_class,
    da_results = da_results,
    metric_name = "Spearman",
    min_abundance = min_abundance
  )

  model_specific_plots <- list(
    shap_spearman = setNames(vector("list", length(shap_models_to_plot)), shap_models_to_plot),
    shap_beeswarm = setNames(vector("list", length(shap_models_to_plot)), shap_models_to_plot),
    shap_prevalence = setNames(vector("list", length(shap_models_to_plot)), shap_models_to_plot),
    biomarker_integrated_spearman = setNames(vector("list", length(shap_models_to_plot)), shap_models_to_plot)
  )

  if (length(shap_models_to_plot) > 0L) {
    for (model_name in shap_models_to_plot) {
      model_specific_plots$shap_spearman[[model_name]] <- plotShapGlobal(
        importance$shap_top$spearman,
        outcome_var,
        class_of_interest = class_of_interest,
        negative_class = negative_class,
        top_n = top_n,
        metric_name = "Spearman",
        model = model_name
      )

      model_specific_plots$shap_beeswarm[[model_name]] <- plotShapBeeswarm(
        importance$shap_raw,
        top_n = top_n,
        model = model_name
      )
      model_specific_plots$shap_prevalence[[model_name]] <- plotTaxaPrevalence(
        filtered_counts = plot_filtered_counts,
        metadata = plot_metadata,
        global_importance = importance$shap_top$spearman,
        target_var = outcome_var,
        min_abundance = min_abundance,
        top_n = top_n,
        model = model_name
      )
      model_specific_plots$biomarker_integrated_spearman[[model_name]] <- plotBiomarkerIntegrated(
        filtered_counts = plot_filtered_counts,
        metadata = plot_metadata,
        global_importance = importance$shap_top$spearman,
        target_var = outcome_var,
        top_n = top_n,
        class_of_interest = class_of_interest,
        negative_class = negative_class,
        da_results = da_results,
        metric_name = "Spearman",
        min_abundance = min_abundance,
        model = model_name
      )

    }
  }

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

    shap_beeswarm = plotShapBeeswarm(
      importance$shap_raw, top_n = top_n
    ),
    shap_prevalence = plotTaxaPrevalence(
      filtered_counts = plot_filtered_counts,
      metadata = plot_metadata,
      global_importance = importance$global_importance$spearman,
      target_var = outcome_var,
      min_abundance = min_abundance,
      top_n = top_n
    ),
    biomarker_integrated_spearman = plotBiomarkerIntegrated(
      filtered_counts = plot_filtered_counts,
      metadata = plot_metadata,
      global_importance = importance$global_importance$spearman,
      target_var = outcome_var,
      top_n = top_n,
      class_of_interest = class_of_interest,
      negative_class = negative_class,
      da_results = da_results,
      metric_name = "Spearman",
      min_abundance = min_abundance
    ),

    perm_heat = plotPermutationHeatmap(
      importance$permutation_top
    ),
    model_specific = model_specific_plots
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
    plot_names <- c("metrics", "roc", "cm", "shap_spearman",
                    "shap_beeswarm", "shap_prevalence",
                    "biomarker_integrated_spearman", "perm_heat")

    # Define sizes (width, height)
    plot_sizes <- list(
      metrics                       = c(8, 5),
      roc                           = c(8, 5),
      cm                            = c(10, 8),
      shap_spearman                 = c(9, 8),
      shap_beeswarm                 = c(10, 8),
      shap_prevalence               = c(9, 8),
      biomarker_integrated_spearman = c(12, 9),
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

    if (length(shap_models_to_plot) > 0L) {
      for (model_name in shap_models_to_plot) {
        model_slug <- gsub("[^A-Za-z0-9]+", "_", model_name)

        ggplot2::ggsave(
          filename = paste0(prefix, "_", model_slug, "_shap_spearman.png"),
          plot = plots$model_specific$shap_spearman[[model_name]],
          width = 9, height = 8, dpi = 300
        )

        ggplot2::ggsave(
          filename = paste0(prefix, "_", model_slug, "_shap_beeswarm.png"),
          plot = plots$model_specific$shap_beeswarm[[model_name]],
          width = 10, height = 8, dpi = 300
        )
        ggplot2::ggsave(
          filename = paste0(prefix, "_", model_slug, "_shap_prevalence.png"),
          plot = plots$model_specific$shap_prevalence[[model_name]],
          width = 9, height = 8, dpi = 300
        )
        ggplot2::ggsave(
          filename = paste0(prefix, "_", model_slug, "_biomarker_integrated_spearman.png"),
          plot = plots$model_specific$biomarker_integrated_spearman[[model_name]],
          width = 12, height = 9, dpi = 300
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
    # >>> ADICIONE ESTAS LINHAS ABAIXO <<<
    readr::write_tsv(
      integrated_table_consensus,
      paste0(prefix, "_integrated_biomarkers.tsv")
    )

    if (run_da && !is.null(da_results)) {
      # Filtra os resultados de DA apenas para os Top Táxons
      top_da_results <- da_results %>%
        dplyr::filter(Taxon %in% integrated_table_consensus$Taxon)

      readr::write_tsv(
        top_da_results,
        paste0(prefix, "_top_differential_abundance.tsv")
      )
    }
    # >>> FIM DA ADIÇÃO <<<

    readr::write_tsv(
      importance$shap_top$spearman,
      paste0(prefix, "_shap_model_spearman.tsv")
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
    models          = final_models,
    metrics         = metrics_summary,
    importance      = importance,
    plots           = plots,
    data            = list(train = train_data, test = test_data),
    selected_models = selected_shap_models,
    da_results      = da_results,
    integrated_table = integrated_table_consensus # <--- Adicione esta linha
  )
}
