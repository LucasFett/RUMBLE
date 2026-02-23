#' Plotting Functions for RUMBLE Results
#'
#' A collection of visualization functions for model performance,
#' feature importance, and interpretability results.
#'
#' @name plot-functions
#' @keywords internal
#' @return No return value. This object is internal and serves only for organization.
#'
#' @import ggplot2
#' @importFrom ggsci scale_fill_nejm scale_color_nejm
#' @importFrom viridis scale_fill_viridis
#' @importFrom rlang sym
#' @importFrom dplyr group_by summarize mutate left_join pull filter
#'      slice_head arrange desc ungroup select
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom purrr imap_dfr
#' @importFrom yardstick roc_curve roc_auc conf_mat
NULL


#' Plot ROC Curves for All Models
#'
#' Generates ROC curves with AUC values in the legend for each model
#' evaluated on the test set.
#'
#' @param final_models Named list of model results from the pipeline.
#' @param target_var Character. Name of the target variable.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'   # Mock data structure
#'   mock_preds <- data.frame(
#'     Status = factor(rep(c("Control", "Disease"), each = 10)),
#'     .pred_Disease = runif(20),
#'     .pred_class = factor(sample(c("Control", "Disease"), 20, replace = TRUE))
#'   )
#'
#'   mock_models <- list(
#'     RF = list(predictions = mock_preds),
#'     XGB = list(predictions = mock_preds)
#'   )
#'
#'   plotRocCurves(mock_models, "Status")
#'
#'
#' @export
plotRocCurves <- function(final_models, target_var) {
  roc_data <- purrr::imap_dfr(final_models, function(obj, name) {
    df <- obj$predictions
    pos_class <- levels(df[[target_var]])[2L]
    prob_col  <- paste0(".pred_", pos_class)

    roc_tbl <- yardstick::roc_curve(
      df,
      truth = !!rlang::sym(target_var),
      !!rlang::sym(prob_col),
      event_level = "second"
    )

    auc_val <- yardstick::roc_auc(
      df,
      truth = !!rlang::sym(target_var),
      !!rlang::sym(prob_col),
      event_level = "second"
    )$.estimate

    roc_tbl$model <- name
    roc_tbl$auc   <- auc_val
    roc_tbl
  })

  auc_labels <- roc_data %>%
    dplyr::group_by(model) %>%
    dplyr::summarize(
      AUC = unique(round(auc, 2L)),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = paste0(model, " (AUC = ", AUC, ")")
    )

  roc_data <- dplyr::left_join(roc_data, auc_labels, by = "model")

  ggplot2::ggplot(
    roc_data,
    ggplot2::aes(
      x = 1 - specificity, y = sensitivity,
      color = label
    )
  ) +
    ggplot2::geom_path(linewidth = 1) +
    ggplot2::geom_abline(
      linetype = "dashed", color = "gray50"
    ) +
    ggsci::scale_color_nejm() +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = "ROC Curves",
      subtitle = paste("Target:", target_var),
      x = "1 - Specificity", y = "Sensitivity",
      color = "Model"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Plot Global SHAP Consensus
#'
#' Produces a directional barplot of the top features ranked by mean
#' absolute SHAP value across all models. Color indicates the direction
#' of the effect.
#'
#' @param global_importance A data.frame from
#'      \code{\link{computeFeatureImportance}} containing columns
#'      \code{variable}, \code{mean_shap}, and \code{direction}.
#' @param target_var Character. Name of the target variable.
#' @param positive_class Character. Label of the positive class.
#' @param negative_class Character. Label of the negative class.
#' @param top_n Integer. Number of top features to display (default 20).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'   # Mock data
#'   mock_global <- data.frame(
#'     variable = paste0("Taxon_", 1:5),
#'     mean_shap = runif(5, 0.1, 0.5),
#'     direction = sample(c(-1, 1), 5, replace = TRUE)
#'   )
#'
#'   plotShapGlobal(mock_global, "Status", "Disease", "Control", top_n = 5)
#'
#'
#' @export
plotShapGlobal <- function(global_importance, target_var,
                           positive_class, negative_class,
                           top_n = 20L) {
  df_plot <- head(global_importance, top_n)

  df_plot$directional_value <- ifelse(
    df_plot$direction < 0,
    df_plot$mean_shap * -1,
    df_plot$mean_shap
  )

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = reorder(variable, mean_shap),
      y = directional_value, fill = direction
    )
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient2(
      low = "#d62728", mid = "#f7f7f7", high = "#1f77b4",
      midpoint = 0, limits = c(-1, 1),
      name = "Effect direction"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Mean SHAP contribution (directional)",
      title = "Biomarker Consensus (SHAP)",
      subtitle = paste0(
        target_var, ": <- ", negative_class,
        " | ", positive_class, " ->"
      )
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Plot Confusion Matrices
#'
#' Displays confusion matrices as heatmaps with absolute counts and
#' proportions for each model.
#'
#' @param final_models Named list of model results.
#' @param target_var Character. Name of the target variable.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotConfusionMatrices <- function(final_models, target_var) {
  cm_df <- purrr::imap_dfr(final_models, function(x, name) {
    cm <- yardstick::conf_mat(
      x$predictions,
      truth = !!rlang::sym(target_var),
      estimate = .pred_class
    )
    df <- as.data.frame(cm$table)
    df$model <- name
    df
  })

  cm_df <- cm_df %>%
    dplyr::group_by(model, Truth) %>%
    dplyr::mutate(Prop = Freq / sum(Freq) * 100) %>%
    # Lógica para contraste: Se a barra for escura (>50%), texto branco
    dplyr::mutate(text_col = ifelse(Prop > 50, "white", "black")) %>%
    dplyr::ungroup()

  ggplot2::ggplot(
    cm_df,
    ggplot2::aes(x = Prediction, y = Truth, fill = Prop)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(Freq, "\n(", round(Prop, 1), "%)"),
        color = text_col  # Usa a cor calculada dinamicamente
      ),
      size = 3.5, fontface = "bold"
    ) +
    ggplot2::facet_wrap(~ model) +
    # Paleta Gradiente Azul (Clean e Profissional)
    ggplot2::scale_fill_gradient(
      low = "#eff3ff",  # Azul muito claro (quase branco)
      high = "#08519c", # Azul escuro forte
      limits = c(0, 100),
      name = "Proportion (%)"
    ) +
    # Garante que o ggplot entenda que "white" e "black" são cores reais
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      title = "Confusion Matrices (Test Set)",
      x = "Predicted", y = "Actual"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black") # Eixos pretos para leitura
    )
}


#' Plot Metrics Comparison Barplot
#'
#' Produces a grouped barplot comparing performance metrics (F1, AUC,
#' Accuracy, Balanced Accuracy) across all models.
#'
#' @param final_models Named list of model results.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'   # Mock metrics
#'   mock_met <- data.frame(
#'     .metric = c("accuracy", "roc_auc"),
#'     .estimate = c(0.85, 0.92)
#'   )
#'   mock_models <- list(RF = list(metrics = mock_met))
#'
#'   plotMetricsComparison(mock_models)
#'
#'
#' @export
plotMetricsComparison <- function(final_models) {
  metrics_df <- purrr::imap_dfr(final_models, function(obj, name) {
    obj$metrics$model <- name
    obj$metrics
  })

  metrics_long <- metrics_df %>%
    dplyr::select(model, .metric, .estimate) %>%
    tidyr::pivot_wider(
      names_from = .metric, values_from = .estimate
    ) %>%
    tidyr::pivot_longer(
      cols = -model, names_to = "Metric",
      values_to = "Value"
    )

  ggplot2::ggplot(
    metrics_long,
    ggplot2::aes(x = model, y = Value, fill = Metric)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.7
    ) +
    ggsci::scale_fill_nejm() +
    ggplot2::labs(
      title = "Model Performance Comparison",
      x = "Model", y = "Metric Value", fill = "Metric"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::ylim(0, 1) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1
      )
    )
}


#' Plot SHAP Distribution by Model
#'
#' Displays a jitter + boxplot showing the distribution of SHAP values
#' for the top features, broken down by model.
#'
#' @param shap_raw A data.frame of raw SHAP values from
#'      \code{\link{computeFeatureImportance}}.
#' @param top_features Character vector of feature names to display.
#'      If \code{NULL}, the top \code{top_n} features are computed
#'      automatically.
#' @param top_n Integer. Number of top features to show if
#'      \code{top_features} is \code{NULL} (default 15).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'   # Mock SHAP raw data
#'   mock_shap <- data.frame(
#'     variable = rep(c("Taxon_A", "Taxon_B"), each = 10),
#'     contribution = rnorm(20),
#'     model = rep(c("RF", "XGB"), 10)
#'   )
#'
#'   plotShapDistribution(mock_shap, top_n = 2)
#'
#'
#' @export
plotShapDistribution <- function(shap_raw, top_features = NULL,
                                 top_n = 15L) {
  if (is.null(top_features)) {
    top_features <- shap_raw %>%
      dplyr::group_by(variable) %>%
      dplyr::summarize(
        mean_abs = mean(abs(contribution)),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(mean_abs)) %>%
      dplyr::slice_head(n = top_n) %>%
      dplyr::pull(variable)
  }

  shap_filtered <- shap_raw[shap_raw$variable %in% top_features, ]

  ggplot2::ggplot(
    shap_filtered,
    ggplot2::aes(
      y = reorder(variable, abs(contribution)),
      x = contribution,
      color = sign(contribution)
    )
  ) +
    ggplot2::geom_jitter(
      size = 1.5, alpha = 0.6,
      ggplot2::aes(shape = model)
    ) +
    ggplot2::geom_boxplot(
      alpha = 0.4, fill = "orange", color = "black",
      outlier.alpha = 0
    ) +
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "black", linewidth = 0.5
    ) +
    ggplot2::labs(
      y = "Feature", x = "SHAP Contribution",
      title = "SHAP Distribution by Model",
      shape = "Model", color = "Direction"
    ) +
    ggplot2::scale_color_gradient2(
      low = "#d62728", mid = "gray", high = "#1f77b4",
      midpoint = 0
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Plot Permutation Importance Barplot
#'
#' Displays a grouped barplot of the top features ranked by
#' permutation-based variable importance (delta loss), per model.
#'
#' @param permutation_top A data.frame of top permutation features from
#'      \code{\link{computeFeatureImportance}}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'   # Mock data
#'   mock_perm <- data.frame(
#'     variable = rep(c("Taxon_A", "Taxon_B"), 2),
#'     mean_dropout_loss = c(0.05, 0.02, 0.04, 0.01),
#'     model = rep(c("RF", "XGB"), each = 2)
#'   )
#'
#'   plotPermutationImportance(mock_perm)
#'
#'
#' @export
plotPermutationImportance <- function(permutation_top) {
  # Filter > 0 to match user logic, but ensure variable exists
  df_plot <- permutation_top[permutation_top$mean_dropout_loss > 0, ]

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      y = reorder(variable, mean_dropout_loss),
      x = mean_dropout_loss, fill = model
    )
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(
      y = "Feature",
      x = "Importance (delta loss)",
      title = "Permutation Importance by Model",
      fill = "Model"
    ) +
    ggsci::scale_fill_nejm() +
    ggplot2::theme_minimal(base_size = 12)
}


#' Plot Permutation Importance Heatmap
#'
#' Displays a heatmap of permutation-based variable importance for each
#' feature across models.
#'
#' @param permutation_top A data.frame of top permutation features from
#'      \code{\link{computeFeatureImportance}}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#'   # Mock data
#'   mock_perm <- data.frame(
#'     variable = rep(c("Taxon_A", "Taxon_B"), 2),
#'     mean_dropout_loss = c(0.05, 0.02, 0.04, 0.01),
#'     model = rep(c("RF", "XGB"), each = 2)
#'   )
#'   plotPermutationHeatmap(mock_perm)
#'
#'
#' @export
plotPermutationHeatmap <- function(permutation_top) {
  df_plot <- permutation_top[permutation_top$mean_dropout_loss > 0, ]

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      y = reorder(variable, mean_dropout_loss),
      x = model, fill = mean_dropout_loss
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::labs(
      y = "Feature", x = "Model",
      fill = "Importance\n(delta loss)",
      title = "Permutation Importance Heatmap"
    ) +
    viridis::scale_fill_viridis(option = "H") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1
      )
    )
}


#' Summarise Model Metrics
#'
#' Extracts and pivots performance metrics from all models into a
#' single summary table.
#'
#' @param final_models Named list of model results.
#'
#' @return A data.frame with one row per model and one column per
#'      metric.
#'
#' @examples
#'   # Mock data
#'   mock_met <- data.frame(
#'     .metric = c("accuracy", "roc_auc"),
#'     .estimate = c(0.85, 0.92)
#'   )
#'   mock_models <- list(RF = list(metrics = mock_met))
#'
#'   summariseMetrics(mock_models)
#'
#'
#' @export
summariseMetrics <- function(final_models) {
  purrr::imap_dfr(final_models, function(obj, name) {
    obj$metrics$model <- name
    obj$metrics
  }) %>%
    dplyr::select(model, .metric, .estimate) %>%
    tidyr::pivot_wider(
      names_from = .metric, values_from = .estimate
    )
}
