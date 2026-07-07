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
#' @importFrom stats reorder setNames
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

#' Select and Standardize Importance Table
#'
#' Internal helper to filter importance tables by model and standardize the
#' importance column name for plotting functions.
#'
#' @param importance_df Data.frame of feature importances.
#' @param model Character or NULL. Model name to filter by.
#' @return Data.frame with standardized `importance_value` column.
#' @noRd
.selectImportanceTable <- function(importance_df, model = NULL) {
  df <- importance_df

  # Filtra pelo modelo, se solicitado
  if (!is.null(model)) {
    if (!"model" %in% colnames(df)) {
      stop("Cannot filter by model: 'model' column not found in importance table.")
    }
    df <- df[df$model == model, , drop = FALSE]
  }

  # Padroniza o nome da coluna de importância para os gráficos
  if ("mean_abs_contribution" %in% colnames(df)) {
    df$importance_value <- df$mean_abs_contribution
  } else if ("contribution" %in% colnames(df)) {
    df$importance_value <- abs(df$contribution)
  } else {
    stop("Could not find a valid importance column (expected 'mean_abs_contribution').")
  }

  # Padroniza a direção (Rho)
  if ("SHAP_Rho" %in% colnames(df)) {
    df$direction <- df$SHAP_Rho
  } else if (!"direction" %in% colnames(df)) {
    df$direction <- 0
  }

  return(df)
}

#' Plot Cross-Validated ROC Curves with Variance Band
#'
#' Generates publication-quality ROC curves evaluating model performance across
#' cross-validation folds. For each model, it computes the empirical ROC curve
#' and AUC for each fold individually, then averages sensitivities across a standard
#' false positive rate grid to display a robust mean ROC curve surrounded by a
#' standard deviation (SD) variance ribbon. Individual fold paths are displayed
#' as thin, semi-transparent background lines.
#'
#' @param final_models Named list of model results from the pipeline (consolidated out-of-fold structure).
#' @param target_var Character. Name of the target variable.
#'
#' @return A \code{ggplot} object faceted by model.
#'
#' @export
plotRocCurves <- function(final_models, target_var) {

  grid_fpr <- seq(0, 1, by = 0.01)

  fold_roc_list <- list()
  auc_summary_list <- list()

  for (model_name in names(final_models)) {
    obj <- final_models[[model_name]]
    df <- obj$predictions

    if (!is.factor(df[[target_var]])) {
      df[[target_var]] <- as.factor(df[[target_var]])
    }

    pos_class <- levels(df[[target_var]])[2L]
    prob_col  <- paste0(".pred_", pos_class)

    folds <- unique(df$fold)
    if (is.null(folds)) folds <- 1

    for (f in folds) {
      df_fold <- df[df$fold == f, , drop = FALSE]

      if (length(unique(df_fold[[target_var]])) < 2 || nrow(df_fold) < 2) next

      # Calcula a curva ROC específica deste fold
      curve_fold <- yardstick::roc_curve(
        df_fold,
        truth = !!rlang::sym(target_var),
        !!rlang::sym(prob_col),
        event_level = "second"
      )

      # Calcula a AUC específica deste fold
      auc_fold <- yardstick::roc_auc(
        df_fold,
        truth = !!rlang::sym(target_var),
        !!rlang::sym(prob_col),
        event_level = "second"
      )$.estimate

      fpr_raw <- 1 - curve_fold$specificity
      tpr_raw <- curve_fold$sensitivity

      # Ordenação rigorosa para garantir o funcionamento da interpolação linear
      ord <- order(fpr_raw)
      fpr_raw <- fpr_raw[ord]
      tpr_raw <- tpr_raw[ord]

      # Interpola as sensibilidades para o grid padronizado de FPR
      interp_tpr <- stats::approx(x = fpr_raw, y = tpr_raw, xout = grid_fpr, ties = max)$y

      fold_roc_list[[length(fold_roc_list) + 1]] <- data.frame(
        model = model_name,
        fold = f,
        fpr = grid_fpr,
        tpr = interp_tpr
      )

      auc_summary_list[[length(auc_summary_list) + 1]] <- data.frame(
        model = model_name,
        fold = f,
        auc = auc_fold
      )
    }
  }

  if (length(fold_roc_list) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No ROC data available"))
  }

  df_folds_roc <- dplyr::bind_rows(fold_roc_list)
  df_auc <- dplyr::bind_rows(auc_summary_list)

  # Agrega os resultados para calcular a média e desvio padrão ponto a ponto
  df_mean_roc <- df_folds_roc %>%
    dplyr::group_by(model, fpr) %>%
    dplyr::summarize(
      mean_tpr = mean(tpr, na.rm = TRUE),
      sd_tpr = sd(tpr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Restringe os limites matemáticos para não estourarem o espaço biológico [0,1]
      ymin = pmax(0, mean_tpr - sd_tpr),
      ymax = pmin(1, mean_tpr + sd_tpr)
    )

  # Gera as métricas de sumário (Média ± SD) da AUC entre os folds externos
  df_auc_stats <- df_auc %>%
    dplyr::group_by(model) %>%
    dplyr::summarize(
      mean_auc = mean(auc, na.rm = TRUE),
      sd_auc = sd(auc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      label = paste0(model, " (AUC = ", sprintf("%.3f", mean_auc), " \u00b1 ", sprintf("%.3f", sd_auc), ")")
    )

  df_mean_roc <- dplyr::left_join(df_mean_roc, df_auc_stats, by = "model")
  df_folds_roc <- dplyr::left_join(df_folds_roc, df_auc_stats, by = "model")

  p <- ggplot2::ggplot() +
    # Linha diagonal de referência (classificador aleatório)
    ggplot2::geom_abline(linetype = "dashed", color = "gray60", linewidth = 0.5) +
    # Banda de variância (Sombreamento do Desvio Padrão)
    ggplot2::geom_ribbon(
      data = df_mean_roc,
      ggplot2::aes(x = fpr, ymin = ymin, ymax = ymax, fill = model),
      alpha = 0.15, inherit.aes = FALSE
    ) +
    # Caminhos individuais de cada fold (linhas finas e claras em background)
    ggplot2::geom_line(
      data = df_folds_roc,
      ggplot2::aes(x = fpr, y = tpr, group = interaction(model, fold), color = model),
      linewidth = 0.4, alpha = 0.35
    ) +
    # Curva média consolidada do modelo (linha espessa em destaque)
    ggplot2::geom_line(
      data = df_mean_roc,
      ggplot2::aes(x = fpr, y = mean_tpr, color = model),
      linewidth = 1.2
    ) +
    ggplot2::facet_wrap(~ label) +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    # Pequena folga (-0.01 a 1.01) para resolver definitivamente o corte visual nas bordas
    ggplot2::scale_x_continuous(limits = c(-0.01, 1.01), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(limits = c(-0.01, 1.01), expand = c(0, 0), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(
      title = "Cross-Validated ROC Curves by Model",
      subtitle = paste("Target:", target_var, "| Shaded area represents mean \u00b1 SD across folds"),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.spacing = ggplot2::unit(1.2, "lines"),
      strip.text = ggplot2::element_text(face = "bold", size = 10.5),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::coord_fixed()

  return(p)
}

#' Plot Global SHAP Consensus or Model-Specific SHAP Profiles
#'
#' Produces a directional barplot of the top features ranked by SHAP
#' importance. By default, it uses consensus rankings across models, but it can
#' also display the isolated profile of a specific model when a per-model SHAP
#' table is supplied.
#'
#' @param global_importance A data.frame from
#'      \code{\link{computeFeatureImportance}} containing columns
#'      \code{variable}, \code{direction}, and \code{mean_abs_contribution}.
#' @param target_var Character. Name of the target variable.
#' @param class_of_interest Character. Label of the class of interest.
#' @param negative_class Character. Label of the negative/reference class.
#' @param top_n Integer. Number of top features to display (default 20).
#' @param metric_name Character. Name of the metric used (e.g., "Spearman").
#' @param model Character or \code{NULL}. Optional model name for isolated
#'   plots when \code{global_importance} contains a \code{model} column.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotShapGlobal <- function(global_importance, target_var,
                           class_of_interest,
                           negative_class,
                           top_n = 20L,
                           metric_name = "Spearman",
                           model = NULL) {
  if (missing(class_of_interest)) {
    stop("'class_of_interest' is required in plotShapGlobal().")
  }

  if (missing(negative_class)) {
    stop("'negative_class' is required in plotShapGlobal().")
  }

  df_plot <- .selectImportanceTable(global_importance, model = model)
  df_plot <- df_plot %>%
    dplyr::arrange(dplyr::desc(importance_value)) %>%
    dplyr::slice_head(n = top_n)

  df_plot$directional_value <- ifelse(
    df_plot$direction < 0,
    df_plot$importance_value * -1,
    df_plot$importance_value
  )

  df_plot$label_with_sig <- df_plot$variable

  plot_title <- if (is.null(model)) {
    paste0("Biomarker Consensus (", metric_name, ")")
  } else {
    paste0("Biomarker Profile: ", model, " (", metric_name, ")")
  }

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = stats::reorder(label_with_sig, importance_value),
      y = directional_value, fill = direction
    )
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient2(
      low = "#1f77b4", mid = "gray85", high = "#d62728",
      midpoint = 0, limits = c(-1, 1),
      name = "Correlation (\u03c1)"
    ) +
    ggplot2::labs(
      x = NULL,
      y = if (is.null(model)) paste0("Normalized SHAP Score (", metric_name, ")") else paste0("Raw SHAP Score (", metric_name, ")"),
      title = plot_title,
      # Inserida quebra de linha (\n) para evitar o corte do texto
      subtitle = paste0(
        target_var, ": <- ", negative_class,
        " | ", class_of_interest, " ->\n(Direction: Spearman \u03c1 between abundance\nand SHAP impact)"
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 10) # Margem de segurança
    )
}

#' Plot Confusion Matrices (Pooled Out-of-Fold)
#'
#' Displays confusion matrices as heatmaps with absolute counts and
#' proportions for each model using pooled Out-of-Fold (OOF) predictions.
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
      title = "Pooled Out-of-Fold Confusion Matrices",
      subtitle = "Aggregated predictions from external validation folds across Nested CV",
      x = "Predicted", y = "Actual"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "black"),
      plot.title = ggplot2::element_text(face = "bold")
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
#' for the top features using ggbeeswarm.
#'
#' @param shap_raw A data.frame of raw SHAP values.
#' @param prediction_data Data.frame used for predictions, containing the target variable. Optional.
#' @param target_var Character. Name of the target variable. Optional.
#' @param top_n Integer. Number of top features to show (default 20).
#' @param model Character or \code{NULL}. Optional model name to isolate a
#'   specific model in the beeswarm plot.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotShapBeeswarm <- function(shap_raw, prediction_data = NULL, target_var = NULL, top_n = 20L, model = NULL) {

  if (!is.null(model)) {
    available_models <- unique(shap_raw$model)
    if (!all(model %in% available_models)) {
      stop(
        "Invalid model selection. Available models: ",
        paste(available_models, collapse = ", ")
      )
    }
    shap_raw <- shap_raw[shap_raw$model %in% model, , drop = FALSE]
  }

  # Obter as top features ranqueadas pela média absoluta do SHAP
  top_features <- shap_raw %>%
    dplyr::group_by(variable) %>%
    dplyr::summarize(mean_abs = mean(abs(contribution), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abs)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(variable)

  shap_filtered <- shap_raw %>%
    dplyr::filter(variable %in% top_features)

  # Normalizar feature values (Min-Max scaling por taxon)
  shap_filtered <- shap_filtered %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(
      feature_val_scaled = (feature_value - min(feature_value, na.rm = TRUE)) /
        (max(feature_value, na.rm = TRUE) - min(feature_value, na.rm = TRUE) + 1e-6)
    ) %>%
    dplyr::ungroup()

  # Fazer o join com as classes reais
  has_classes <- !is.null(prediction_data) && !is.null(target_var)
  if (has_classes) {
    class_df <- data.frame(
      observation_id = seq_len(nrow(prediction_data)),
      True_Class = prediction_data[[target_var]]
    )
    shap_filtered <- dplyr::left_join(shap_filtered, class_df, by = "observation_id")
  }

  plot_title <- if (is.null(model)) {
    "SHAP Distribution by Model"
  } else {
    paste0("SHAP Distribution: ", model)
  }

  # Plot base
  p <- ggplot2::ggplot(
    shap_filtered,
    ggplot2::aes(
      x = contribution,
      y = reorder(variable, abs(contribution), FUN = mean)
    )
  ) +
    # Linha zero central
    ggplot2::geom_vline(xintercept = 0, color = "black")

  # Adicionar os pontos com algoritmo swarm e corral gutter
  if (has_classes) {
    p <- p + ggbeeswarm::geom_beeswarm(
      ggplot2::aes(color = feature_val_scaled, shape = True_Class),
      size = 1.5,
      method = "swarm",
      corral = "gutter",
      corral.width = 0.7,
      priority = "density",
      cex = 0.5
    ) +
      ggplot2::labs(shape = "Class")
  } else {
    p <- p + ggbeeswarm::geom_beeswarm(
      ggplot2::aes(color = feature_val_scaled),
      size = 1.5,
      method = "swarm",
      corral = "gutter",
      corral.width = 0.7,
      priority = "density",
      cex = 0.5
    )
  }

  # Aplicar o gradiente Viridis e o tema clássico
  p <- p + viridis::scale_color_viridis(
    option = "H",
    alpha = 0.7,
    breaks = c(0, 1),
    labels = c("Low", "High"),
    name = "Feature Value"
  ) +
    ggplot2::labs(
      x = "SHAP Value (Impact on Model Output)",
      y = "Feature",
      title = plot_title
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "italic")
    )

  # Se for o plot de consenso, separar em painéis
  if (is.null(model) && length(unique(shap_filtered$model)) > 1) {
    p <- p + ggplot2::facet_wrap(~ model, scales = "free_x")
  }

  return(p)
}

#' Plot Taxa Prevalence for Top SHAP Features
#'
#' Displays the prevalence of top SHAP features in the dataset.
#' This function uses raw/relative abundance data (filtered_counts) to calculate
#' true prevalence, avoiding artifacts from CLR transformation.
#'
#' @param filtered_counts A data.frame of raw or relative abundances (before CLR/rCLR
#'   transformation) with samples as rows and taxa as columns. This should be the
#'   \code{filtered_counts} element from \code{prepareData()}.
#' @param metadata A data.frame with sample metadata, including the target variable.
#'   Rownames should match the sample identifiers in \code{filtered_counts}.
#' @param global_importance Global importance data.frame or per-model importance
#'   table.
#' @param target_var Character. Name of the target variable.
#' @param min_abundance Numeric. Minimum relative abundance threshold to
#'   consider a taxon as present in a sample (default 0.0001). Not used directly if
#'   filtered_counts are raw counts, but kept for compatibility.
#' @param top_n Integer. Number of top features to show (default 20).
#' @param model Character or \code{NULL}. Optional model name for isolated
#'   prevalence plots when \code{global_importance} contains a \code{model}
#'   column.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotTaxaPrevalence <- function(filtered_counts, metadata, global_importance,
                               target_var, min_abundance = 0.0001,
                               top_n = 20L, model = NULL) {
  top_taxa_df <- .selectImportanceTable(global_importance, model = model) %>%
    dplyr::arrange(dplyr::desc(importance_value)) %>%
    dplyr::slice_head(n = top_n)

  top_taxa <- top_taxa_df$variable
  taxa_order <- top_taxa_df$variable[order(top_taxa_df$importance_value)]

  # Ensure filtered_counts has the top taxa
  available_taxa <- intersect(top_taxa, colnames(filtered_counts))
  if (length(available_taxa) == 0) {
    stop("No top taxa found in filtered_counts. Check taxon names.")
  }

  # Extract subset of data
  X_data <- filtered_counts[, available_taxa, drop = FALSE]
  groups <- metadata[[target_var]]

  # CRITICAL FIX: Ensure groups is a factor so levels() works
  if (!is.factor(groups)) {
    groups <- factor(groups)
  }

  # Ensure alignment
  if (length(groups) != nrow(X_data)) {
    stop("Number of samples in filtered_counts and metadata must match.")
  }

  # Calculate prevalence
  prev_list <- lapply(levels(groups), function(lvl) {
    idx <- which(groups == lvl)
    X_group <- X_data[idx, , drop = FALSE]

    presence_matrix <- X_group > 0
    prev <- colMeans(presence_matrix)

    data.frame(variable = names(prev),
               prevalence = as.numeric(prev),
               group = lvl,
               stringsAsFactors = FALSE)
  })

  df_prev <- dplyr::bind_rows(prev_list)
  df_prev$variable <- factor(df_prev$variable, levels = taxa_order)
  class_colors <- stats::setNames(c("#1f77b4", "#d62728"), levels(groups))

  plot_title <- if (is.null(model)) {
    "Prevalence of Top Biomarkers"
  } else {
    paste0("Prevalence of Top Biomarkers: ", model)
  }

  ggplot2::ggplot(
    df_prev,
    ggplot2::aes(x = prevalence, y = variable, fill = group)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::scale_fill_manual(values = class_colors) +
    ggplot2::labs(
      title = plot_title,
      x = "Prevalence (Fraction of Samples)",
      y = NULL,
      fill = target_var
    ) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    ggplot2::theme_minimal(base_size = 12)
}


#' Plot Integrated Biomarker Interpretation (Dashboard Style)
#'
#' Generates a 4-panel visualization dissecting the consensus biomarkers:
#' \itemize{
#'   \item \strong{Panel A:} ML Importance (Median Abs SHAP, colored by direction).
#'   \item \strong{Panel B:} Taxa Prevalence across groups.
#'   \item \strong{Panel C:} Biological Effect Size (Log2 Fold-Change).
#'   \item \strong{Panel D:} Statistical Significance (-log10 FDR).
#' }
#'
#' @param filtered_counts A numeric matrix of counts (taxa as columns).
#' @param metadata A data.frame with sample metadata.
#' @param global_importance A data.frame with consensual feature importance.
#' @param target_var Character. Name of the outcome variable.
#' @param top_n Integer. Number of top features to include.
#' @param class_of_interest Character. Label of the class of interest.
#' @param negative_class Character. Label of the reference class.
#' @param da_results A data.frame with differential abundance results. Optional.
#' @param metric_name Character. Feature importance metric name (e.g., "Spearman").
#' @param min_abundance Numeric. Threshold for prevalence calculation.
#' @param model Character. Optional model name to plot model-specific SHAP.
#' @param verbose Logical. Print messages.
#'
#' @return A patchworked ggplot object.
#' @export
plotBiomarkerIntegrated <- function(filtered_counts,
                                    metadata,
                                    global_importance,
                                    target_var,
                                    top_n = 20,
                                    class_of_interest,
                                    negative_class,
                                    da_results = NULL,
                                    metric_name = "Spearman",
                                    min_abundance = 0.0001,
                                    model = NULL,
                                    verbose = TRUE) {

  msg <- function(...) { if (verbose) message(...) }
  msg("Generating integrated biomarker dashboard plot...")

  show_da <- !is.null(da_results)
  method_full <- ifelse(is.null(da_results$Method[1]), "DA", da_results$Method[1])

  # ==================================================================
  # 1. Preparar dados
  # ==================================================================
  importance_full <- .selectImportanceTable(global_importance, model = model)

  if (!"direction" %in% colnames(importance_full)) importance_full$direction <- 0

  importance_full <- importance_full %>%
    dplyr::mutate(label_with_sig = variable)

  taxa_data_top <- importance_full %>%
    dplyr::filter(variable != "TOTAL_SHAP") %>%
    dplyr::slice_max(order_by = importance_value, n = top_n, with_ties = FALSE) %>%
    dplyr::arrange(importance_value)

  taxa_data_top$label_with_sig <- factor(taxa_data_top$label_with_sig, levels = taxa_data_top$label_with_sig)

  taxa_order_labels <- levels(taxa_data_top$label_with_sig)
  taxa_order_original <- taxa_data_top$variable

  sig_map <- stats::setNames(taxa_data_top$label_with_sig, taxa_data_top$variable)

  ## ==================================================================
  ## Painel A: ML Importance
  ## ==================================================================
  msg("  - Generating Panel A (ML Importance with directional bars)...")

  taxa_data_top <- taxa_data_top %>%
    dplyr::mutate(
      directional_value = ifelse(direction < 0, importance_value * -1, importance_value)
    )

  max_imp <- max(abs(taxa_data_top$directional_value), na.rm = TRUE)

  p_imp <- ggplot2::ggplot(taxa_data_top, ggplot2::aes(x = directional_value, y = label_with_sig, fill = direction)) +
    ggplot2::geom_col(color = "black", linewidth = 0.3, width = 0.7) +
    ggplot2::scale_fill_gradient2(
      low = "#1f77b4", mid = "gray85", high = "#d62728",
      midpoint = 0, limits = c(-1, 1),
      name = paste0(metric_name, " (\u03c1)")
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    ggplot2::labs(
      title = "A. ML Importance",
      subtitle = if (is.null(model)) "Normalized Directional SHAP" else "Raw Directional SHAP",
      x = "Importance",
      y = NULL
    ) +
    ggplot2::scale_x_continuous(limits = c(-max_imp * 1.05, max_imp * 1.05)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9, face = "italic"),
                   panel.grid.major.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9, color = "gray30"))

  ## ==================================================================
  ## Painel B: Taxa Prevalence
  ## ==================================================================
  msg("  - Generating Panel B (Taxa Prevalence)...")

  groups <- metadata[[target_var]]
  rel_abund <- sweep(filtered_counts, 1, rowSums(filtered_counts) + 1e-9, "/")

  avail_taxa <- intersect(taxa_order_original, colnames(rel_abund))

  df_prev <- lapply(levels(factor(groups)), function(lvl) {
    idx <- which(groups == lvl)
    prev <- colMeans(rel_abund[idx, avail_taxa, drop = FALSE] > min_abundance)
    data.frame(variable = names(prev), prevalence = as.numeric(prev), group = lvl)
  }) %>% dplyr::bind_rows()

  df_prev$label_with_sig <- sig_map[df_prev$variable]
  df_prev$label_with_sig <- factor(df_prev$label_with_sig, levels = taxa_order_labels)

  p_prev <- ggplot2::ggplot(df_prev, ggplot2::aes(x = prevalence, y = label_with_sig, fill = group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c("#1f77b4", "#d62728"), name = "Group") +
    ggplot2::scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
    ggplot2::labs(title = "B. Prevalence", subtitle = paste0("Abund > ", scales::percent(min_abundance)), x = "Fraction", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 11, face = "bold"))

  ## ==================================================================
  ## Painéis C & D: Differential Abundance (Reintegrado)
  ## ==================================================================
  if (show_da) {
    da_plot_data <- da_results %>%
      dplyr::filter(Taxon %in% taxa_order_original) %>%
      dplyr::mutate(
        label_with_sig = factor(sig_map[Taxon], levels = taxa_order_labels),
        neg_log_fdr = -log10(ifelse(adj_p_val == 0, 1e-16, adj_p_val)),
        fc_direction = ifelse(logFC > 0, "Positive", "Negative")
      )

    p_fc <- ggplot2::ggplot(da_plot_data, ggplot2::aes(x = logFC, y = label_with_sig, fill = fc_direction)) +
      ggplot2::geom_col(color = "black", linewidth = 0.3, width = 0.7) +
      ggplot2::scale_fill_manual(values = c("Positive" = "#d62728", "Negative" = "#1f77b4"), guide = "none") +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::labs(title = "C. Effect Size", subtitle = method_full, x = "Log2FC", y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = 11, face = "bold"))

    p_sig <- ggplot2::ggplot(da_plot_data, ggplot2::aes(x = neg_log_fdr, y = label_with_sig)) +
      ggplot2::geom_col(fill = "gray80", color = "black", linewidth = 0.3, width = 0.7) +
      ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
      ggplot2::labs(title = "D. Significance", subtitle = "-log10 FDR", x = "Score", y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = 11, face = "bold"))

    plot_list <- list(p_imp, p_prev, p_fc, p_sig)
    plot_layout <- patchwork::plot_layout(widths = c(1.5, 1, 1, 0.8), guides = "collect")
  } else {
    p_fc <- ggplot2::ggplot() + ggplot2::theme_void()
    plot_list <- list(p_imp, p_prev, p_fc)
    plot_layout <- patchwork::plot_layout(widths = c(1.5, 1, 1), guides = "collect")
  }

  ## ==================================================================
  ## Assemblage
  ## ==================================================================
  msg("Assembling dashboard plot...")

  annotation_title <- if (is.null(model)) {
    "Integrated Biomarker Dashboard"
  } else {
    paste0("Integrated Biomarker Dashboard: ", model)
  }

  composite_plot <- patchwork::wrap_plots(plot_list) + plot_layout +
    patchwork::plot_annotation(
      title = annotation_title,
      subtitle = "Direction: Spearman correlation (\u03c1) between abundance and SHAP impact",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold", margin = ggplot2::margin(b = 10)))
    ) &
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "horizontal",
                   legend.margin = ggplot2::margin(t = 10, b = 10),
                   legend.text = ggplot2::element_text(size = 11),
                   legend.title = ggplot2::element_text(size = 12, face = "bold"))

  return(composite_plot)
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

#' Extract Integrated Biomarker Data
#'
#' Builds a comprehensive data.frame containing SHAP importance, taxonomic
#' prevalence, and differential abundance metrics for the top features. This
#' table mirrors the data shown in the Integrated Biomarker Dashboard plot.
#'
#' @inheritParams plotBiomarkerIntegrated
#'
#' @return A data.frame with the integrated metrics.
#' @export
getIntegratedBiomarkerTable <- function(filtered_counts,
                                        metadata,
                                        global_importance,
                                        target_var,
                                        top_n = 20,
                                        class_of_interest,
                                        negative_class,
                                        da_results = NULL,
                                        metric_name = "Spearman",
                                        min_abundance = 0.0001,
                                        model = NULL) {

  # 1. SHAP Importance
  importance_full <- .selectImportanceTable(global_importance, model = model)

  taxa_order <- importance_full %>%
    dplyr::filter(variable != "TOTAL_SHAP") %>%
    dplyr::slice_max(order_by = importance_value, n = top_n, with_ties = FALSE) %>%
    dplyr::arrange(dplyr::desc(importance_value)) %>%
    dplyr::pull(variable)

  df_integrated <- importance_full %>%
    dplyr::filter(variable %in% taxa_order) %>%
    dplyr::mutate(
      Taxon = variable,
      Direction_Class = ifelse(direction > 0, class_of_interest, negative_class),
      SHAP_Score = importance_value,
      Direction = direction
    ) %>%
    dplyr::select(Taxon, SHAP_Score, Direction_Class, Direction)

  # 2. Prevalence
  groups <- metadata[[target_var]]
  depths <- rowSums(filtered_counts)
  depths[depths == 0] <- 1
  rel_abund <- sweep(filtered_counts, 1, depths, "/")

  avail_taxa <- intersect(df_integrated$Taxon, colnames(rel_abund))
  X_data <- rel_abund[, avail_taxa, drop = FALSE]

  prev_list <- lapply(levels(factor(groups)), function(lvl) {
    idx <- which(groups == lvl)
    X_group <- X_data[idx, , drop = FALSE]
    presence_matrix <- X_group > min_abundance
    prev <- colMeans(presence_matrix)
    data.frame(Taxon = names(prev),
               Prevalence = as.numeric(prev),
               Group = lvl,
               stringsAsFactors = FALSE)
  })

  df_prev <- dplyr::bind_rows(prev_list) %>%
    tidyr::pivot_wider(names_from = Group, values_from = Prevalence, names_prefix = "Prevalence_")

  # Merge SHAP and Prevalence
  df_integrated <- dplyr::left_join(df_integrated, df_prev, by = "Taxon")

  # 3. Differential Abundance (If available)
  if (!is.null(da_results)) {
    da_subset <- da_results %>%
      dplyr::filter(Taxon %in% df_integrated$Taxon) %>%
      dplyr::select(Taxon, logFC, p_val, adj_p_val) %>%
      dplyr::rename(DA_logFC = logFC, DA_p_val = p_val, DA_adj_p_val = adj_p_val)

    df_integrated <- dplyr::left_join(df_integrated, da_subset, by = "Taxon")
  } else {
    df_integrated$DA_logFC <- NA
    df_integrated$DA_p_val <- NA
    df_integrated$DA_adj_p_val <- NA
  }

  # Ensure order matches the plot (descending importance)
  df_integrated <- df_integrated %>%
    dplyr::arrange(match(Taxon, taxa_order))

  return(df_integrated)
}



#' Plot Out-of-Fold Predicted Probability Density
#'
#' Displays the distribution of predicted probabilities from the held-out
#' test folds, separated by true class. This reveals calibration issues and
#' class separability that aggregate metrics cannot capture.
#'
#' @param final_models Named list of model results with predictions.
#' @param target_var Character. Name of the target variable.
#'
#' @return A \code{ggplot} object faceted by model.
#'
#' @export
plotOOFDensity <- function(final_models, target_var) {
  # Collect all OOF predictions
  df_all <- purrr::imap_dfr(final_models, function(obj, name) {
    preds <- obj$predictions
    pos_class <- levels(preds[[target_var]])[2L]
    prob_col <- paste0(".pred_", pos_class)

    if (!prob_col %in% colnames(preds)) return(NULL)

    data.frame(
      model = name,
      probability = preds[[prob_col]],
      true_class = as.character(preds[[target_var]]),
      stringsAsFactors = FALSE
    )
  })

  if (nrow(df_all) == 0) return(NULL)

  ggplot2::ggplot(df_all, ggplot2::aes(x = probability, fill = true_class)) +
    ggplot2::geom_density(alpha = 0.5, color = "black", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dashed", color = "gray40") +
    ggplot2::facet_wrap(~ model, scales = "free_y") +
    ggplot2::scale_fill_manual(values = c("#1f77b4", "#d62728")) +
    ggplot2::labs(
      title = "Out-of-Fold Predicted Probability Density",
      subtitle = "Separation between classes indicates discriminative power. Dashed line = 0.5 threshold.",
      x = "Predicted Probability (Class of Interest)",
      y = "Density",
      fill = "True Class"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold", size = 11)
    )
}


#' Plot SHAP vs Permutation Importance Scatter
#'
#' Visualizes the relationship between SHAP importance (credit distribution)
#' and Permutation Importance (global loss). Features that rank high in both
#' methods are the most robust biomarkers.
#'
#' @param global_importance Data.frame of global SHAP importance.
#' @param permutation_top Data.frame of top permutation features.
#' @param feature_frequency Data.frame of feature selection frequency (optional).
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotSHAPvsPermutation <- function(global_importance, permutation_top, feature_frequency = NULL) {

  # Agrega a permutação (tirando a média entre os modelos)
  perm_agg <- permutation_top %>%
    dplyr::group_by(variable) %>%
    dplyr::summarize(
      perm_importance = mean(mean_delta_loss, na.rm = TRUE),
      .groups = "drop"
    )

  # Pega o SHAP global
  shap_df <- global_importance %>%
    dplyr::select(variable, mean_abs_contribution) %>%
    dplyr::rename(shap_importance = mean_abs_contribution)

  # Junta as duas métricas
  df_merged <- dplyr::inner_join(shap_df, perm_agg, by = "variable")

  if (nrow(df_merged) == 0) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No overlapping features for SHAP vs Permutation"))

  # Define as medianas para traçar os quadrantes
  med_shap <- stats::median(df_merged$shap_importance, na.rm = TRUE)
  med_perm <- stats::median(df_merged$perm_importance, na.rm = TRUE)

  # Marca quem está no quadrante superior direito (os mais robustos)
  df_merged <- df_merged %>%
    dplyr::mutate(
      is_top_quadrant = (shap_importance > med_shap) & (perm_importance > med_perm),
      point_color = ifelse(is_top_quadrant, "#d62728", "gray60"),
      label = ifelse(is_top_quadrant, variable, "") # Só coloca nome nos top
    )

  p <- ggplot2::ggplot(df_merged, ggplot2::aes(x = perm_importance, y = shap_importance)) +
    # Linhas dos quadrantes
    ggplot2::geom_vline(xintercept = med_perm, linetype = "dashed", color = "gray80") +
    ggplot2::geom_hline(yintercept = med_shap, linetype = "dashed", color = "gray80") +
    # Pontos fixos (destacando em vermelho os top)
    ggplot2::geom_point(ggplot2::aes(color = point_color), size = 3, alpha = 0.8) +
    ggplot2::scale_color_identity() +
    # Linha de tendência
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "gray30", linetype = "dashed", linewidth = 0.5) +
    ggplot2::labs(
      title = "SHAP vs Permutation Importance",
      subtitle = "Highlighting robust biomarkers (Upper-Right Quadrant)",
      x = "Permutation Importance (Mean Delta Loss)",
      y = "SHAP Importance (Mean |SHAP|)"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  # Aplica o ggrepel se estiver instalado (o que evita sobreposição de textos)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 3.5,
      fontface = "italic",
      box.padding = 0.5,
      max.overlaps = Inf,
      color = "black"
    )
  } else {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label),
      size = 3.5,
      hjust = -0.1,
      vjust = -0.3,
      check_overlap = TRUE
    )
  }

  return(p)
}
