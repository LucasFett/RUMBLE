# Declare global variables to prevent R CMD check notes when using tidyverse NSE (Non-Standard Evaluation)
utils::globalVariables(c(
  ".estimate", ".metric", ".pred_class", "AUC", "Direction", "Freq",
  "Metric", "Prediction", "Prop", "Taxon", "Truth", "Value", "auc",
  "contribution", "delta_loss", "direction", "directional_value",
  "dropout_loss", "feature_value", "full_model_loss", "group",
  "importance_value", "label", "logFC", "mean_abs", "mean_abs_contribution",
  "mean_delta_loss", "model", "observation_id", "prevalence",
  "sensitivity", "significance", "specificity", "text_col", "variable",
  "Direction_Class", "SHAP_Score", "Prevalence", "Group", "p_val", "adj_p_val",
  "direction_class", "plot_shap_value", "direction_discrete", "neg_log_fdr","Parameter", "Value", "mean", "std_err", "CV_Mean", "CV_StdErr"
))
