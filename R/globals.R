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
  "sd_abs_contribution"
))
