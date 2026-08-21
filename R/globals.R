#' @importFrom stats sd
NULL

# Declare global variables to prevent R CMD check notes when using tidyverse NSE (Non-Standard Evaluation)
utils::globalVariables(c(
  # yardstick / tidymodels
  ".estimate", ".metric", ".pred_class", "AUC",

  # Plot aesthetics and columns
  "Direction", "Direction_Class", "Freq", "Metric", "Prediction", "Prop",
  "Taxon", "Truth", "Value", "auc", "directional_value", "group",
  "importance_value", "label", "label_with_sig", "logFC", "mean_abs",
  "mean_abs_contribution", "mean_delta_loss", "model", "observation_id",
  "prevalence", "sensitivity", "specificity", "text_col", "variable",
  "SHAP_Score", "Prevalence", "Group",
  "direction_class", "plot_shap_value", "direction_discrete", "neg_log_fdr",
  "fc_direction",

  # explain_functions.R
  "contribution", "delta_loss", "direction", "dropout_loss", "feature_value",
  "full_model_loss", "weight", "imp",

  # run_analysis.R
  "Parameter", "CV_Mean", "CV_StdErr", "Mean", "SD", "Fold",
  "n_folds_present", "fold_frequency",

  # DA columns
  "p_val", "adj_p_val", "DA_logFC", "DA_p_val", "DA_adj_p_val",

  # plotOOFDensity
  "probability", "true_class",

  # plotSHAPvsPermutation
  "perm_importance", "shap_importance",

  # hyperparameter stability
  "n_unique_values", "stability_flag",

  # ROC
  "fpr", "tpr_mean", "tpr_min", "tpr_max", "mean_auc",

  # feature frequency
  "sd_abs_contribution",

  # --- Adições da Auditoria Final (Variáveis Auxiliares e Métricas) ---
  ".config", "Accuracy_Rate", "Best_CV_Mean", "Discrepancy_Max_Min",
  "Formatted_Value", "Generalization_Gap", "Is_Correct", "Kendall",
  "Magnitude Change (%)", "Magnitude_Change_Pct", "Mean_Inner", "Mean_Outer",
  "Mean_SHAP_Test", "Mean_SHAP_Train", "Model", "Models_Failed",
  "Prob_True_Class", "SD_Inner", "Spearman", "Spearman (Top N)",
  "True_Class", "abs_contribution", "feature_val_scaled", "fold",
  "importance", "is_top_quadrant", "mean_abs_shap_test", "mean_abs_shap_train",
  "is_stable",
  "mean_estimate", "mean_tpr", "model_freq", "model_n_folds", "point_color",
  "sample_id", "scaled_abs_contribution", "sd_auc", "sd_estimate", "sd_tpr",
  "std_err", "tpr", "ymax", "ymin"
))
