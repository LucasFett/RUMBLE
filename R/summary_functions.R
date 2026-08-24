#' Calculate Feature Selection Frequency Across Folds (Consensus and Per-Model)
#'
#' Calculates the fraction of times a feature appeared in the top-N features.
#' Treats each model-fold combination as an independent evaluation. For example,
#' 4 models across 4 folds equals 16 total evaluations.
#'
#' @param shap_raw Data.frame of raw SHAP values with 'fold' and 'model' columns.
#' @param top_n Integer. Number of top features per fold.
#' @param n_folds Integer. Total number of outer folds.
#' @return Data.frame with democratic global and per-model frequencies.
#' @export
calculateFeatureFrequency <- function(shap_raw, top_n = 20L, n_folds = 5L, min_fold_frequency = 0.7) {
  if (nrow(shap_raw) == 0 || !"fold" %in% colnames(shap_raw)) {
    return(data.frame(variable = character(0), n_folds_present = integer(0), fold_frequency = numeric(0)))
  }

  models_present <- unique(shap_raw$model)
  n_models <- length(models_present)
  total_evaluations <- as.numeric(n_models * n_folds)

  # 1. Extracts the Top N individually for EACH model within EACH fold
  model_fold_tops <- shap_raw %>%
    dplyr::group_by(model, fold, variable) %>%
    dplyr::summarize(imp = mean(abs(contribution), na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(model, fold) %>%
    dplyr::slice_max(order_by = imp, n = top_n, with_ties = FALSE) %>%
    dplyr::ungroup()

  # 2. Global Frequency (Democratic Count)
  freq_global <- model_fold_tops %>%
    dplyr::group_by(variable) %>%
    dplyr::summarize(
      n_folds_present = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      total_possible_evaluations = total_evaluations,
      fold_frequency = as.numeric(n_folds_present) / total_evaluations
    )

  # 3. Per-Model Specific Frequency
  freq_model_wide <- model_fold_tops %>%
    dplyr::group_by(model, variable) %>%
    dplyr::summarize(
      model_n_folds = dplyr::n_distinct(fold),
      .groups = "drop"
    ) %>%
    dplyr::mutate(model_freq = as.numeric(model_n_folds) / as.numeric(n_folds)) %>%
    tidyr::pivot_wider(
      names_from = model,
      values_from = c(model_n_folds, model_freq),
      values_fill = 0,
      names_glue = "{model}_{.value}"
    )

  # 4. Final Join (Consensus + Per Model)
  result <- freq_global %>%
    dplyr::full_join(freq_model_wide, by = "variable") %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~tidyr::replace_na(., 0))) %>%
    dplyr::mutate(is_stable = fold_frequency >= min_fold_frequency) %>%
    dplyr::arrange(dplyr::desc(fold_frequency), dplyr::desc(n_folds_present))

  return(result)
}
#' Summarize Outer CV Performance Consistency
#'
#' Evaluates the stability of model performance across the outer folds
#' of the Nested Cross-Validation, identifying potential high-variance models.
#'
#' @param metrics_per_fold Data.frame containing metrics per model per fold.
#' @param metrics_to_report Character vector of metrics to include in the formatted table.
#' @return A list containing the full raw stats and the publication-ready formatted table.
#' @export
summarizeOuterCV <- function(metrics_per_fold, metrics_to_report = c("roc_auc", "mcc", "bal_accuracy")) {

  if (is.null(metrics_per_fold) || nrow(metrics_per_fold) == 0) return(NULL)

  # 1. Base Calculations
  stats_df <- metrics_per_fold %>%
    dplyr::select(model, .metric, fold, .estimate) %>%
    tidyr::pivot_wider(names_from = fold, values_from = .estimate, names_prefix = "Fold_") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Mean = mean(dplyr::c_across(dplyr::starts_with("Fold_")), na.rm = TRUE),
      SD = sd(dplyr::c_across(dplyr::starts_with("Fold_")), na.rm = TRUE),
      Discrepancy_Max_Min = max(dplyr::c_across(dplyr::starts_with("Fold_")), na.rm = TRUE) -
        min(dplyr::c_across(dplyr::starts_with("Fold_")), na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.metric, model)

  # 2. Formatting for Publication
  final_table <- stats_df %>%
    dplyr::filter(.metric %in% metrics_to_report) %>%
    dplyr::mutate(
      Formatted_Value = sprintf("%.3f \u00B1 %.3f", Mean, SD),
      Discrepancy_Max_Min = round(Discrepancy_Max_Min, 3), # Keeps it as NUMERIC
      Model = dplyr::case_when(
        model == "RF" ~ "Random Forest (RF)",
        model == "XGB" ~ "XGBoost (XGB)",
        model == "ENET" ~ "Elastic Net (ENET)",
        model == "KNN" ~ "KNN",
        TRUE ~ model
      )
    ) %>%
    dplyr::select(Model, .metric, Formatted_Value, Discrepancy_Max_Min) %>%
    tidyr::pivot_wider(
      names_from = .metric,
      values_from = c(Formatted_Value, Discrepancy_Max_Min),
      names_glue = "{.name} ({.value})"
    )

  # 3. Column Name Cleanup
  colnames_to_fix <- colnames(final_table)
  colnames_to_fix <- gsub("Formatted_Value", "Outer CV Mean \u00B1 SD", colnames_to_fix)
  colnames_to_fix <- gsub("Discrepancy_Max_Min", "Max Discrepancy", colnames_to_fix)
  colnames_to_fix <- gsub("_", " ", colnames_to_fix)
  colnames(final_table) <- toupper(colnames_to_fix)
  colnames(final_table)[1] <- "Model"

  # Ordering
  auc_col <- grep("ROC AUC \\(OUTER CV MEAN", colnames(final_table), value = TRUE)
  if (length(auc_col) > 0) {
    final_table <- final_table %>% dplyr::arrange(dplyr::desc(!!rlang::sym(auc_col[1])))
  }

  return(list(raw_stats = stats_df, publication_table = final_table))
}

#' Summarize Inner CV Tuning Consistency
#'
#' Evaluates the stability of the winning hyperparameters across the inner folds.
#'
#' @param cv_metrics Data.frame containing inner CV grid search results.
#' @param metrics_to_report Character vector of metrics to include.
#' @return A list containing the full raw stats and the publication-ready formatted table.
#' @export
summarizeInnerCV <- function(cv_metrics, metrics_to_report = c("roc_auc", "mcc")) {

  if (is.null(cv_metrics) || nrow(cv_metrics) == 0) return(NULL)

  # Since cv_metrics now arrives already filtered at the source
  # (run_analysis.R), we can select the values directly without applying max()
  best_inner <- cv_metrics %>%
    dplyr::select(Model, Fold, Metric, Best_CV_Mean = CV_Mean) %>%
    dplyr::distinct() # Extra safeguard against accidental duplicate rows

  stats_inner <- best_inner %>%
    tidyr::pivot_wider(names_from = Fold, values_from = Best_CV_Mean, names_prefix = "Inner_Fold_") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      Mean_Inner = mean(dplyr::c_across(dplyr::starts_with("Inner_Fold_")), na.rm = TRUE),
      SD_Inner = sd(dplyr::c_across(dplyr::starts_with("Inner_Fold_")), na.rm = TRUE),
      Discrepancy_Max_Min = max(dplyr::c_across(dplyr::starts_with("Inner_Fold_")), na.rm = TRUE) -
        min(dplyr::c_across(dplyr::starts_with("Inner_Fold_")), na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Metric, Model)

  final_table_inner <- stats_inner %>%
    dplyr::filter(Metric %in% metrics_to_report) %>%
    dplyr::mutate(
      Formatted_Value = sprintf("%.3f \u00B1 %.3f", Mean_Inner, SD_Inner),
      Discrepancy_Max_Min = round(Discrepancy_Max_Min, 3), # Keeps it as NUMERIC
      Model = dplyr::case_when(
        Model == "RF" ~ "Random Forest (RF)",
        Model == "XGB" ~ "XGBoost (XGB)",
        Model == "ENET" ~ "Elastic Net (ENET)",
        Model == "KNN" ~ "KNN",
        TRUE ~ Model
      )
    ) %>%
    dplyr::select(Model, Metric, Formatted_Value, Discrepancy_Max_Min) %>%
    tidyr::pivot_wider(
      names_from = Metric,
      values_from = c(Formatted_Value, Discrepancy_Max_Min),
      names_glue = "{.name} ({.value})"
    )

  colnames_to_fix <- colnames(final_table_inner)
  colnames_to_fix <- gsub("Formatted_Value", "Inner CV Mean \u00B1 SD", colnames_to_fix)
  colnames_to_fix <- gsub("Discrepancy_Max_Min", "Max Discrepancy", colnames_to_fix)
  colnames_to_fix <- gsub("_", " ", colnames_to_fix)
  colnames(final_table_inner) <- toupper(colnames_to_fix)
  colnames(final_table_inner)[1] <- "Model"

  auc_col <- grep("ROC AUC \\(INNER CV MEAN", colnames(final_table_inner), value = TRUE)
  if (length(auc_col) > 0) {
    final_table_inner <- final_table_inner %>% dplyr::arrange(dplyr::desc(!!rlang::sym(auc_col[1])))
  }

  return(list(raw_stats = stats_inner, publication_table = final_table_inner))
}

#' Analyze Out-of-Fold Predictions (Confidence & Hard Samples)
#'
#' Evaluates the prediction confidence (probability assigned to the true class)
#' for correct vs. incorrect predictions, and identifies "hard samples" that
#' were consistently misclassified across all surviving models.
#'
#' @param oof_preds Data.frame of out-of-fold predictions.
#' @param target_var Character. Name of the target variable.
#' @return A list containing confidence summary and hard samples data.frames.
#' @export
analyzeOOFPredictions <- function(oof_preds, target_var) {
  if (is.null(oof_preds) || nrow(oof_preds) == 0) return(NULL)

  # Ensures the truth column is text, to make matching against the .pred_ columns easier
  truth_vec <- as.character(oof_preds[[target_var]])

  # 1. Mean Confidence Analysis
  # Dynamically extracts the probability assigned to the TRUE class of that sample
  probs <- vapply(seq_len(nrow(oof_preds)), function(i) {
    class_name <- truth_vec[i]
    prob_col <- paste0(".pred_", class_name)

    if (prob_col %in% colnames(oof_preds)) {
      return(oof_preds[[prob_col]][i])
    } else {
      return(NA_real_)
    }
  }, numeric(1L)) # <-- This is the key part: we declare the return type as numeric

  conf_df <- oof_preds %>%
    dplyr::mutate(
      Truth = truth_vec,
      Is_Correct = as.character(.pred_class) == Truth,
      Prob_True_Class = probs
    )

  confidence_summary <- conf_df %>%
    dplyr::group_by(model, Is_Correct) %>%
    dplyr::summarize(
      N_Samples = dplyr::n(),
      Mean_Prob_True_Class = mean(Prob_True_Class, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(model, dplyr::desc(Is_Correct))

  # 2. Identification of Hard Samples
  # Evaluates samples against all models present in the OOF preds
  hard_samples <- conf_df %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarize(
      True_Class = dplyr::first(Truth),
      Accuracy_Rate = mean(Is_Correct, na.rm = TRUE),
      Models_Failed = paste(model[!Is_Correct], collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::filter(Accuracy_Rate == 0) %>% # Failed in 100% of the available models
    dplyr::select(sample_id, True_Class, Models_Failed) %>%
    dplyr::arrange(True_Class)

  return(list(
    confidence = confidence_summary,
    hard_samples = hard_samples
  ))
}
#' Evaluate SHAP Magnitude Stability (Train vs Test Drop)
#'
#' Calculates the percentage change in SHAP magnitude between the training
#' and test sets for the Top-N features (defined by training importance).
#' This measures "Explainability Overfitting" - a collapse in magnitude
#' indicates the model's learned biological rules do not generalize well.
#'
#' @param shap_raw_train Data.frame of raw SHAP values from the training set.
#' @param shap_raw_test Data.frame of raw SHAP values from the test set.
#' @param top_n Integer. Number of top features from the train set to evaluate.
#' @return A formatted data.frame summarizing the magnitude drop per model.
#' @export
evaluateShapMagnitudeStability <- function(shap_raw_train, shap_raw_test, top_n = 20L) {
  if (is.null(shap_raw_train) || nrow(shap_raw_train) == 0 ||
      is.null(shap_raw_test)  || nrow(shap_raw_test) == 0) {
    return(NULL)
  }

  # 1. Compute the mean importance (absolute SHAP) per feature/model in both phases
  imp_train <- shap_raw_train %>%
    dplyr::group_by(model, variable) %>%
    dplyr::summarize(mean_abs_shap_train = mean(abs(contribution), na.rm = TRUE), .groups = "drop")

  imp_test <- shap_raw_test %>%
    dplyr::group_by(model, variable) %>%
    dplyr::summarize(mean_abs_shap_test = mean(abs(contribution), na.rm = TRUE), .groups = "drop")

  # 2. Join the data to compute the change
  magnitude_drop <- imp_train %>%
    dplyr::inner_join(imp_test, by = c("model", "variable"))

  # 3. Global summary focused only on the Top N features from TRAINING
  resumo_magnitude <- magnitude_drop %>%
    dplyr::group_by(model) %>%
    dplyr::slice_max(order_by = mean_abs_shap_train, n = top_n, with_ties = FALSE) %>%
    dplyr::summarize(
      Mean_SHAP_Train = mean(mean_abs_shap_train, na.rm = TRUE),
      Mean_SHAP_Test = mean(mean_abs_shap_test, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Percentage change: (Final - Initial) / Initial * 100
      Magnitude_Change_Pct = ((Mean_SHAP_Test - Mean_SHAP_Train) / Mean_SHAP_Train) * 100
    ) %>%
    # Using round() to keep a pure NUMERIC format (avoids bugs in spreadsheet software)
    dplyr::mutate(
      Mean_SHAP_Train = round(Mean_SHAP_Train, 5),
      Mean_SHAP_Test = round(Mean_SHAP_Test, 5),
      Magnitude_Change_Pct = round(Magnitude_Change_Pct, 2),
      Model = dplyr::case_when(
        model == "RF" ~ "Random Forest (RF)",
        model == "XGB" ~ "XGBoost (XGB)",
        model == "ENET" ~ "Elastic Net (ENET)",
        model == "KNN" ~ "KNN",
        TRUE ~ model
      )
    ) %>%
    # Arranging columns for the final table
    dplyr::select(
      Model,
      `Mean SHAP (Train Top N)` = Mean_SHAP_Train,
      `Mean SHAP (Test)` = Mean_SHAP_Test,
      `Magnitude Change (%)` = Magnitude_Change_Pct
    ) %>%
    # Sorts to place the models that best preserved the signal (or had the smallest drop) at the top
    dplyr::arrange(dplyr::desc(`Magnitude Change (%)`))

  return(resumo_magnitude)
}
#' Extract Core Consensus Features
#'
#' Identifies the intersecting top features across all surviving models.
#' A feature must be present in the Top-N of EVERY evaluated model to be
#' considered part of the "Core Consensus".
#'
#' @param shap_top A data.frame of top SHAP features per model (from RUMBLE$importance$shap_top$spearman).
#' @param top_n Integer. Number of top features to consider per model.
#' @return A list containing the core consensus taxa data.frame.
#' @export
extractCoreConsensus <- function(shap_top, top_n = 20L) {
  if (is.null(shap_top) || nrow(shap_top) == 0 || !"model" %in% colnames(shap_top)) {
    return(NULL)
  }

  # 1. Ensures a clean extraction of the Top N features per model
  top_features <- shap_top %>%
    dplyr::group_by(model) %>%
    # Uses mean_abs_contribution, the standardized column generated by the package
    dplyr::slice_max(order_by = mean_abs_contribution, n = top_n, with_ties = FALSE) %>%
    dplyr::ungroup()

  # 2. Splits the features into a list where each element is one model's Top N
  taxa_per_model <- split(top_features$variable, top_features$model)

  # If only 1 model survived the cutoff, the consensus is simply its own Top N
  if (length(taxa_per_model) == 1) {
    return(list(
      core_consensus = data.frame(Taxon = taxa_per_model[[1]], stringsAsFactors = FALSE),
      n_models = 1
    ))
  }

  # 3. Computes the strict intersection (Core Consensus) across ALL models
  core_taxa <- Reduce(intersect, taxa_per_model)

  # 4. Returns a clean data.frame
  return(list(
    core_consensus = data.frame(Taxon = core_taxa, stringsAsFactors = FALSE),
    n_models = length(taxa_per_model),
    models_included = paste(names(taxa_per_model), collapse = ", ")
  ))
}
#' Evaluate Generalization Gap (Predictive Overfitting)
#'
#' Calculates the discrepancy between the Inner CV (Validation) performance
#' and the Outer CV (Test) performance. A large positive gap indicates that
#' the model memorized the training data (Predictive Overfitting) and failed
#' to generalize to unseen data.
#'
#' @param outer_raw Data.frame of raw outer CV stats (from summarizeOuterCV$raw_stats).
#' @param inner_raw Data.frame of raw inner CV stats (from summarizeInnerCV$raw_stats).
#' @param metrics_to_report Character vector of metrics to include.
#' @return A publication-ready data.frame summarizing the generalization gap.
#' @export
evaluateGeneralizationGap <- function(outer_raw, inner_raw, metrics_to_report = c("roc_auc", "mcc")) {
  if (is.null(outer_raw) || nrow(outer_raw) == 0 ||
      is.null(inner_raw) || nrow(inner_raw) == 0) {
    return(NULL)
  }

  # 1. Standardize and extract the Outer CV (Test) means
  outer_clean <- outer_raw %>%
    dplyr::select(Model = model, Metric = .metric, Mean_Outer = Mean)

  # 2. Standardize and extract the Inner CV (Train/Validation) means
  inner_clean <- inner_raw %>%
    dplyr::select(Model, Metric, Mean_Inner)

  # 3. Join and compute the Gap (Inner - Outer)
  gap_df <- inner_clean %>%
    dplyr::inner_join(outer_clean, by = c("Model", "Metric")) %>%
    dplyr::filter(Metric %in% metrics_to_report) %>%
    dplyr::mutate(
      # Delta > 0 means performance was lost on the test set (Overfitting)
      # Delta < 0 means performance was better on the test set (Underfitting or Sample Variance)
      Generalization_Gap = Mean_Inner - Mean_Outer,

      # Ensure a purely numeric format with 3 decimal places
      Mean_Inner = round(Mean_Inner, 3),
      Mean_Outer = round(Mean_Outer, 3),
      Generalization_Gap = round(Generalization_Gap, 3),

      Model = dplyr::case_when(
        Model == "RF" ~ "Random Forest (RF)",
        Model == "XGB" ~ "XGBoost (XGB)",
        Model == "ENET" ~ "Elastic Net (ENET)",
        Model == "KNN" ~ "KNN",
        TRUE ~ Model
      )
    ) %>%
    dplyr::arrange(Metric, dplyr::desc(Generalization_Gap))

  # 4. Pivot into the Publication Table
  final_table <- gap_df %>%
    tidyr::pivot_wider(
      names_from = Metric,
      values_from = c(Mean_Inner, Mean_Outer, Generalization_Gap),
      names_glue = "{.name} ({.value})"
    )

  # 5. Column Name Cleanup
  colnames_to_fix <- colnames(final_table)
  colnames_to_fix <- gsub("Mean_Inner", "Inner CV Mean", colnames_to_fix)
  colnames_to_fix <- gsub("Mean_Outer", "Outer CV Mean", colnames_to_fix)
  colnames_to_fix <- gsub("Generalization_Gap", "Gap (Inner - Outer)", colnames_to_fix)
  colnames_to_fix <- gsub("_", " ", colnames_to_fix)
  colnames(final_table) <- toupper(colnames_to_fix)
  colnames(final_table)[1] <- "Model"

  # Sort in descending order by the ROC AUC Gap (show the worst overfitting at the top)
  auc_col <- grep("ROC AUC \\(GAP", colnames(final_table), value = TRUE)
  if (length(auc_col) > 0) {
    final_table <- final_table %>% dplyr::arrange(dplyr::desc(!!rlang::sym(auc_col[1])))
  }

  return(final_table)
}

#' Evaluate SHAP Rank Stability (Explainability Overfitting)
#'
#' Calculates Spearman and Kendall rank correlations between Train and Test
#' SHAP importance, strictly focused on the Top-N features learned during training.
#' This avoids the "Zero-Inflation Bias" where models artificially agree on
#' hundreds of irrelevant features with zero SHAP value.
#'
#' @param shap_raw_train Data.frame of raw SHAP values from the training set.
#' @param shap_raw_test Data.frame of raw SHAP values from the test set.
#' @param top_n Integer. Number of top features from the train set to evaluate.
#' @return A formatted data.frame summarizing the rank correlations per model.
#' @export
evaluateShapRankStability <- function(shap_raw_train, shap_raw_test, top_n = 20L) {
  if (is.null(shap_raw_train) || nrow(shap_raw_train) == 0 ||
      is.null(shap_raw_test)  || nrow(shap_raw_test) == 0) {
    return(NULL)
  }

  # 1. Compute mean importance on the training set
  imp_train <- shap_raw_train %>%
    dplyr::group_by(model, variable) %>%
    dplyr::summarize(mean_abs_shap_train = mean(abs(contribution), na.rm = TRUE), .groups = "drop")

  # 2. Compute mean importance on the test set
  imp_test <- shap_raw_test %>%
    dplyr::group_by(model, variable) %>%
    dplyr::summarize(mean_abs_shap_test = mean(abs(contribution), na.rm = TRUE), .groups = "drop")

  # 3. Filter the Top N based on training, join with test, and compute correlations
  rank_stability <- imp_train %>%
    dplyr::group_by(model) %>%
    # Isolates ONLY the actual Top N taxa the model judged to be crucial
    dplyr::slice_max(order_by = mean_abs_shap_train, n = top_n, with_ties = FALSE) %>%
    dplyr::inner_join(imp_test, by = c("model", "variable")) %>%
    dplyr::summarize(
      Spearman = cor(mean_abs_shap_train, mean_abs_shap_test, method = "spearman", use = "complete.obs"),
      Kendall = cor(mean_abs_shap_train, mean_abs_shap_test, method = "kendall", use = "complete.obs"),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Keeping it purely numeric
      Spearman = round(Spearman, 3),
      Kendall = round(Kendall, 3),
      Model = dplyr::case_when(
        model == "RF" ~ "Random Forest (RF)",
        model == "XGB" ~ "XGBoost (XGB)",
        model == "ENET" ~ "Elastic Net (ENET)",
        model == "KNN" ~ "KNN",
        TRUE ~ model
      )
    ) %>%
    dplyr::select(Model, `Spearman (Top N)` = Spearman, `Kendall (Top N)` = Kendall) %>%
    dplyr::arrange(dplyr::desc(`Spearman (Top N)`))

  return(rank_stability)
}
