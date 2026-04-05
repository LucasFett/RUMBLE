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
#'      slice_head arrange desc ungroup select n_distinct bind_rows
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom purrr imap_dfr
#' @importFrom yardstick roc_curve roc_auc conf_mat
#' @importFrom ggbeeswarm geom_quasirandom
#' @importFrom patchwork wrap_plots
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
#' @param class_of_interest Character. Label of the class of interest.
#' @param negative_class Character. Label of the negative/reference class.
#' @param top_n Integer. Number of top features to display (default 20).
#' @param metric_name Character. Name of the metric used (e.g., "Spearman" or "Mean SHAP").
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotShapGlobal <- function(global_importance, target_var,
                           class_of_interest,
                           negative_class,
                           top_n = 20L,
                           metric_name = "Spearman") {
  # Validate parameters
  if (missing(class_of_interest)) {
    stop("'class_of_interest' is required in plotShapGlobal().")
  }

  if (missing(negative_class)) {
    stop("'negative_class' is required in plotShapGlobal().")
  }

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
      low = "#1f77b4", mid = "gray", high = "#d62728",
      midpoint = 0, limits = c(-1, 1),
      name = "Effect direction"
    ) +
    ggplot2::labs(
      x = NULL,
      y = paste0("Mean SHAP contribution (", metric_name, ")"),
      title = paste0("Biomarker Consensus (", metric_name, ")"),
      subtitle = paste0(
        target_var, ": <- ", negative_class,
        " | ", class_of_interest, " ->"
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
        color = text_col
      ),
      size = 3.5, fontface = "bold"
    ) +
    ggplot2::facet_wrap(~ model) +
    ggplot2::scale_fill_gradient(
      low = "#eff3ff",
      high = "#08519c",
      limits = c(0, 100),
      name = "Proportion (%)"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      title = "Confusion Matrices (Test Set)",
      x = "Predicted", y = "Actual"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black")
    )
}


#' Plot Metrics Comparison Barplot
#'
#' Produces a grouped barplot comparing performance metrics across all models.
#'
#' @param final_models Named list of model results.
#'
#' @return A \code{ggplot} object.
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


#' Plot SHAP Beeswarm Plot
#'
#' Displays a beeswarm plot showing the distribution of SHAP values
#' for the top features, broken down by model type.
#'
#' @param shap_raw A data.frame of raw SHAP values.
#' @param top_n Integer. Number of top features to show (default 20).
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotShapBeeswarm <- function(shap_raw, top_n = 20L) {
  top_features <- shap_raw %>%
    dplyr::group_by(variable) %>%
    dplyr::summarize(mean_abs = mean(abs(contribution)), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abs)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(variable)

  shap_filtered <- shap_raw[shap_raw$variable %in% top_features, ]

  # Define the exact direction (-1 for negative, 1 for positive) for the color gradient.
  shap_filtered$Direction <- sign(shap_filtered$contribution)

  ggplot2::ggplot(
    shap_filtered,
    ggplot2::aes(
      y = reorder(variable, abs(contribution)),
      x = contribution
    )
  ) +
    # Zero line
    ggplot2::geom_vline(
      xintercept = 0, linetype = "dashed",
      color = "gray50", linewidth = 0.5
    ) +
    # Beeswarm limpo, sem boxplot de fundo, com parâmetros calibrados contra sobreposição
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(shape = model, color = Direction),
      groupOnX = FALSE,
      size = 1.0,        # Reduzido de 1.8 para evitar que os pontos se esmaguem
      alpha = 0.7,       # Transparência levemente aumentada
      width = 0.35,      # Espalhamento vertical aumentado de 0.25 para 0.35
      stroke = 0.2       # Contorno mais fino nos shapes
    ) +
    ggplot2::scale_color_gradient2(
      low = "#1f77b4", mid = "gray", high = "#d62728",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Direction"
    ) +
    ggplot2::labs(
      y = NULL, # Removido o título "Feature" do eixo Y para manter o gráfico mais limpo
      x = "SHAP Contribution",
      title = "SHAP Distribution by Model",
      shape = "Model"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "gray90", linetype = "dotted"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )
}

#' Plot Taxa Prevalence for Top SHAP Features
#'
#' Displays the prevalence of top SHAP features in the dataset.
#'
#' @param data The data.frame used for SHAP (test or train).
#' @param global_importance Global importance data.frame.
#' @param target_var Character. Name of the target variable.
#' @param top_n Integer. Number of top features to show (default 20).
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotTaxaPrevalence <- function(data, global_importance, target_var, top_n = 20L) {
  # 1. Pega os top táxons e cria uma ordem fixa baseada na importância (mean_shap)
  top_taxa_df <- head(global_importance, top_n)
  top_taxa <- top_taxa_df$variable
  taxa_order <- top_taxa_df$variable[order(top_taxa_df$mean_shap)]

  X_data <- data[, !colnames(data) %in% target_var, drop = FALSE]
  groups <- data[[target_var]]

  # CRITICAL FIX: Detect if data is rCLR (has exact zeros) or CLR (no exact zeros)
  # rCLR preserves true zeros as exactly 0, while CLR transforms all values
  # Tolerância para resíduos de ponto flutuante
  has_exact_zeros <- (sum(abs(X_data) < 1e-12) / length(X_data)) > 0.05

  prev_list <- lapply(levels(groups), function(lvl) {
    idx <- which(groups == lvl)
    # Pega TODAS as bactérias da amostra (necessário para a lógica de CLR)
    X_group <- X_data[idx, , drop = FALSE]
    # Pega apenas as bactérias de interesse para plotar
    subset_data <- X_group[, top_taxa, drop = FALSE]

    if (has_exact_zeros) {
      # Logic for rCLR: Valores minúsculos são verdadeiras ausências
      presence_matrix <- abs(subset_data) > 1e-12
    } else {
      # CRITICAL FIX: O "zero" verdadeiro é a abundância mínima da linha INTEIRA,
      full_row_mins <- apply(X_group, 1, min)
      presence_matrix <- sweep(subset_data, 1, full_row_mins + 1e-6, ">")
    }

    prev <- colMeans(presence_matrix)
    data.frame(variable = names(prev), prevalence = prev, group = lvl)
  })

  df_prev <- dplyr::bind_rows(prev_list)

  # 2. Trava a variável como um fator na ordem exata do SHAP
  df_prev$variable <- factor(df_prev$variable, levels = taxa_order)

  # 3. Força as cores para baterem com o SHAP: Azul (referência) e Vermelho (interesse)
  class_colors <- stats::setNames(c("#1f77b4", "#d62728"), levels(groups))

  ggplot2::ggplot(
    df_prev,
    # Removemos o reorder() e usamos diretamente a variável transformada em fator
    ggplot2::aes(x = prevalence, y = variable, fill = group)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_fill_manual(values = class_colors) + # Aplica as cores manuais
    ggplot2::labs(
      title = "Prevalence of Top Biomarkers",
      x = "Prevalence (Fraction of Samples)",
      y = NULL,
      fill = target_var
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Plot Integrated Biomarker Interpretation
#'
#' Combines SHAP importance and prevalence information into a single figure.
#'
#' @param data The data.frame used for SHAP.
#' @param global_importance Global importance data.frame.
#' @param target_var Character. Name of the target variable.
#' @param top_n Integer. Number of top features to show (default 20).
#' @param class_of_interest Character.
#' @param negative_class Character.
#' @param metric_name Character. Name of the metric used (e.g., "Spearman" or "Mean SHAP").
#'
#' @return A \code{patchwork} object.
#'
#' @export
plotBiomarkerIntegrated <- function(data, global_importance, target_var,
                                    top_n = 20L, class_of_interest, negative_class,
                                    metric_name = "Spearman") {

  p1 <- plotShapGlobal(global_importance, target_var, class_of_interest, negative_class, top_n, metric_name) +
    ggplot2::labs(title = paste0("SHAP Consensus (", metric_name, ")"), subtitle = NULL)

  p2 <- plotTaxaPrevalence(data, global_importance, target_var, top_n) +
    ggplot2::labs(title = "Group Prevalence", subtitle = NULL) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank())

  # Combine using patchwork
  combined <- p1 + p2 +
    patchwork::plot_layout(widths = c(1, 0.8), guides = "collect") +
    patchwork::plot_annotation(
      title = "Integrated Biomarker Interpretation",
      subtitle = paste0("Top ", top_n, " features by SHAP importance")
    )

  return(combined)
}


#' Plot Permutation Importance Heatmap
#'
#' Displays a heatmap of permutation-based variable importance for each
#' feature across models.
#'
#' @param permutation_top A data.frame of top permutation features.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotPermutationHeatmap <- function(permutation_top) {
  df_plot <- permutation_top[permutation_top$mean_delta_loss > 0, ]

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      y = reorder(variable, mean_delta_loss),
      x = model, fill = mean_delta_loss
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::labs(
      y = "Feature", x = "Model",
      fill = "Importance\n(delta loss)",
      title = "Permutation Importance Heatmap (Delta Loss)"
    ) +
    viridis::scale_fill_viridis(option = "H") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1
      )
    )
}


#' Plot S-R-I Stacked Decomposition
#'
#' Displays the ecological decomposition of the top biomarkers into
#' independence, synergy and redundancy components.
#'
#' @param sri_profile A data.frame created by the TreeSHAP ecology stage.
#' @param top_n Integer. Number of top taxa to display.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotSriDecomposition <- function(sri_profile, top_n = 20L) {
  if (is.null(sri_profile) || nrow(sri_profile) == 0L) {
    return(NULL)
  }

  df_plot <- utils::head(sri_profile, top_n)
  taxa_levels <- rev(df_plot$variable)

  df_long <- df_plot %>%
    dplyr::select(
      variable,
      independence_value,
      synergy_value,
      redundancy_value,
      ecological_role
    ) %>%
    tidyr::pivot_longer(
      cols = c(independence_value, synergy_value, redundancy_value),
      names_to = "component",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      variable = factor(variable, levels = taxa_levels),
      component = factor(
        component,
        levels = c("independence_value", "synergy_value", "redundancy_value"),
        labels = c("Independence", "Synergy", "Redundancy")
      )
    )

  ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = value, y = variable, fill = component)
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_fill_manual(
      values = c(
        Independence = "#1f77b4",
        Synergy = "#d62728",
        Redundancy = "#2ca02c"
      )
    ) +
    ggplot2::labs(
      title = "TreeSHAP S-R-I Decomposition",
      subtitle = "Global SHAP importance partitioned into ecological components",
      x = "Global Importance (absolute SHAP)",
      y = NULL,
      fill = "Component"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Plot S-R-I Ecological Profile Summary
#'
#' Produces a compact point view of the ecological role and component balance
#' for the top taxa recovered by the TreeSHAP ecology stage.
#'
#' @param sri_profile A data.frame created by the TreeSHAP ecology stage.
#' @param top_n Integer. Number of taxa to display.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotSriEcologicalProfile <- function(sri_profile, top_n = 20L) {
  if (is.null(sri_profile) || nrow(sri_profile) == 0L) {
    return(NULL)
  }

  df_plot <- utils::head(sri_profile, top_n)
  df_plot$variable <- factor(df_plot$variable, levels = rev(df_plot$variable))

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = independence_frac,
      y = variable,
      size = total_importance,
      color = ecological_role
    )
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggplot2::scale_color_manual(values = c(Pathogen = "#d62728", Protector = "#1f77b4")) +
    ggplot2::labs(
      title = "Ecological Profile of Top Biomarkers",
      subtitle = "Node size encodes total importance and color encodes disease association",
      x = "Independence fraction",
      y = NULL,
      color = "Role",
      size = "Total importance"
    ) +
    ggplot2::theme_minimal(base_size = 12)
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
