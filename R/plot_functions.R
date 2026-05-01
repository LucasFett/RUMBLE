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


.selectImportanceTable <- function(importance_df, model = NULL) {
  if (!all(c("variable", "direction") %in% colnames(importance_df))) {
    stop("Importance table must contain at least 'variable' and 'direction' columns.")
  }

  score_col <- NULL
  if ("mean_abs_contribution" %in% colnames(importance_df)) {
    score_col <- "mean_abs_contribution"
  }
  if (is.null(score_col)) {
    stop("Importance table must contain 'mean_abs_contribution'.")
  }

  has_model <- "model" %in% colnames(importance_df)

  if (!has_model) {
    if (!is.null(model) && !identical(model, "consensus")) {
      stop("This importance table does not contain per-model results. Use the consensus table or provide a table with a 'model' column.")
    }
    out <- importance_df
    out$model <- "Consensus"
  } else {
    out <- importance_df
    if (!is.null(model)) {
      available_models <- unique(out$model)
      if (!all(model %in% available_models)) {
        stop(
          "Invalid model selection. Available models: ",
          paste(available_models, collapse = ", ")
        )
      }
      out <- out[out$model %in% model, , drop = FALSE]
    }
  }

  out$importance_value <- out[[score_col]]
  out
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

  plot_title <- if (is.null(model)) {
    paste0("Biomarker Consensus (", metric_name, ")")
  } else {
    paste0("Biomarker Profile: ", model, " (", metric_name, ")")
  }

  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = reorder(variable, importance_value),
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
      y = paste0("SHAP contribution (", metric_name, ")"),
      title = plot_title,
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
#' @param model Character or \code{NULL}. Optional model name to isolate a
#'   specific model in the beeswarm plot.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotShapBeeswarm <- function(shap_raw, top_n = 20L, model = NULL) {

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

  top_features <- shap_raw %>%
    dplyr::group_by(variable) %>%
    dplyr::summarize(mean_abs = mean(abs(contribution)), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abs)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(variable)

  shap_filtered <- shap_raw[shap_raw$variable %in% top_features, ]

  # Define the exact direction (-1 for negative, 1 for positive) for the color gradient.
  shap_filtered$Direction <- sign(shap_filtered$contribution)

  plot_title <- if (is.null(model)) {
    "SHAP Distribution by Model"
  } else {
    paste0("SHAP Distribution: ", model)
  }

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
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(shape = model, color = Direction),
      groupOnX = FALSE,
      size = 1.0,
      alpha = 0.7,
      width = 0.35,
      stroke = 0.2
    ) +
    ggplot2::scale_color_gradient2(
      low = "#1f77b4", mid = "gray", high = "#d62728",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Direction"
    ) +
    ggplot2::labs(
      y = NULL,
      x = "SHAP Contribution",
      title = plot_title,
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
  # 1. Prepare shared taxonomic order and STRICT shared levels
  # ==================================================================
  importance_full <- .selectImportanceTable(global_importance, model = model)

  taxa_order <- importance_full %>%
    dplyr::filter(variable != "TOTAL_SHAP") %>%
    dplyr::slice_max(order_by = importance_value, n = top_n, with_ties = FALSE) %>%
    dplyr::arrange(importance_value) %>%
    dplyr::pull(variable)

  # Força a mesma ordem de níveis para todos os painéis (Garante 1 única legenda)
  shared_levels <- c(class_of_interest, negative_class)
  class_colors_discrete <- stats::setNames(c("#d62728", "#1f77b4"), shared_levels)

  importance_plot_data <- importance_full %>%
    dplyr::filter(variable %in% taxa_order) %>%
    dplyr::mutate(
      variable = factor(variable, levels = taxa_order),
      direction_class = ifelse(direction > 0, class_of_interest, negative_class),
      direction_class = factor(direction_class, levels = shared_levels),
      plot_shap_value = ifelse(direction < 0, importance_value * -1, importance_value)
    )

  ## ==================================================================
  ## Panel A: ML Importance (SHAP) - Solid Bars
  ## ==================================================================
  msg("  - Generating Panel A (ML Importance)...")

  p_imp <- ggplot2::ggplot(importance_plot_data, ggplot2::aes(x = plot_shap_value, y = variable, fill = direction_class)) +
    ggplot2::geom_col(color = "black", linewidth = 0.3, width = 0.7) +
    ggplot2::labs(title = "A. ML Importance", subtitle = paste0("Abs SHAP (", metric_name, ")"), x = "SHAP Score", y = NULL) +
    ggplot2::scale_fill_manual(values = class_colors_discrete, name = target_var, na.translate = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9, face = "italic"),
                   panel.grid.major.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9, color = "gray30"),
                   plot.margin = ggplot2::margin(r = 5))

  ## ==================================================================
  ## Panel B: Taxa Prevalence
  ## ==================================================================
  msg("  - Generating Panel B (Taxa Prevalence)...")

  groups <- metadata[[target_var]]

  depths <- rowSums(filtered_counts)
  depths[depths == 0] <- 1
  rel_abund <- sweep(filtered_counts, 1, depths, "/")

  avail_taxa <- intersect(taxa_order, colnames(rel_abund))
  X_data <- rel_abund[, avail_taxa, drop = FALSE]

  prev_list <- lapply(levels(factor(groups)), function(lvl) {
    idx <- which(groups == lvl)
    X_group <- X_data[idx, , drop = FALSE]
    presence_matrix <- X_group > min_abundance
    prev <- colMeans(presence_matrix)
    data.frame(variable = names(prev),
               prevalence = as.numeric(prev),
               group = lvl,
               stringsAsFactors = FALSE)
  })

  df_prev <- dplyr::bind_rows(prev_list)
  df_prev$variable <- factor(df_prev$variable, levels = taxa_order)

  # Força os níveis do Painel B a serem idênticos aos do Painel A
  df_prev$group <- factor(df_prev$group, levels = shared_levels)

  p_prev <- ggplot2::ggplot(df_prev, ggplot2::aes(x = prevalence, y = variable, fill = group)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1), breaks = c(0, 0.5, 1), expand = c(0,0)) +
    ggplot2::scale_fill_manual(values = class_colors_discrete, name = target_var, na.translate = FALSE) +
    # Aqui alteramos a legenda para deixar explícito que é a abundância relativa
    ggplot2::labs(title = "B. Prevalence", subtitle = paste0("Rel. Abund. > ", scales::percent(min_abundance, accuracy = 0.01)), x = "Fraction") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9, color = "gray30"),
                   plot.margin = ggplot2::margin(r = 5, l = 5))

  ## ==================================================================
  ## Panel C & D: Differential Abundance (LogFC and Significance)
  ## ==================================================================
  if (show_da) {
    msg("  - Generating Panels C & D (Differential Abundance)...")

    da_plot_data <- da_results %>%
      dplyr::filter(Taxon %in% taxa_order) %>%
      dplyr::mutate(Taxon = factor(Taxon, levels = taxa_order))

    da_plot_data <- da_plot_data %>%
      dplyr::mutate(
        direction_discrete = ifelse(logFC > 0, class_of_interest, negative_class),
        # Força os níveis do Painel C a serem idênticos aos do Painel A e B
        direction_discrete = factor(direction_discrete, levels = shared_levels),
        neg_log_fdr = -log10(ifelse(adj_p_val == 0, 1e-16, adj_p_val))
      )

    # Panel C: Effect Size (Log2FC)
    p_fc <- ggplot2::ggplot(da_plot_data, ggplot2::aes(x = logFC, y = Taxon, fill = direction_discrete)) +
      ggplot2::geom_col(color = "black", orientation = "y", linewidth = 0.3, width = 0.7) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      ggplot2::labs(title = "C. Effect Size", subtitle = paste0("Log2FC (", method_full, ")"), x = "Log2 Fold-Change") +
      ggplot2::scale_fill_manual(values = class_colors_discrete, name = target_var, na.translate = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = 11, face = "bold"),
                     plot.subtitle = ggplot2::element_text(size = 9, color = "gray30"),
                     plot.margin = ggplot2::margin(r = 5, l = 5))

    # Panel D: Significance (-log10 FDR)
    p_sig <- ggplot2::ggplot(da_plot_data, ggplot2::aes(x = neg_log_fdr, y = Taxon)) +
      ggplot2::geom_col(fill = "gray80", color = "black", linewidth = 0.3, width = 0.7) +
      ggplot2::geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red") +
      ggplot2::labs(title = "D. Significance", subtitle = "-log10(FDR)", x = "-log10(FDR)") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.1))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     panel.grid.major.y = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = 11, face = "bold"),
                     plot.subtitle = ggplot2::element_text(size = 9, color = "gray30"),
                     plot.margin = ggplot2::margin(l = 5))

  } else {
    p_fc <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "C. Effect Size\n(N/A)")
    p_sig <- ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "D. Significance\n(N/A)")
  }

  ## ==================================================================
  ## Assemblage with Unified Theme
  ## ==================================================================
  msg("Assembling 4-panel dashboard plot...")

  plot_list <- list(p_imp, p_prev, p_fc, p_sig)

  plot_layout <- patchwork::plot_layout(widths = c(1.5, 1, 1, 0.8), guides = "collect")

  annotation_title <- if (is.null(model)) {
    "Integrated Biomarker Dashboard"
  } else {
    paste0("Integrated Biomarker Dashboard: ", model)
  }

  composite_plot <- patchwork::wrap_plots(plot_list) + plot_layout +
    patchwork::plot_annotation(
      title = annotation_title,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 15, face = "bold", margin = ggplot2::margin(b = 10)))
    ) &
    ggplot2::theme(legend.position = "bottom",
                   legend.box = "horizontal",
                   legend.margin = ggplot2::margin(t = 10, b = 10),
                   legend.text = ggplot2::element_text(size = 11),
                   legend.title = ggplot2::element_text(size = 12, face = "bold"))

  return(composite_plot)
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
      SHAP_Score = importance_value
    ) %>%
    dplyr::select(Taxon, SHAP_Score, Direction_Class)

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
      dplyr::select(Taxon, logFC, p_val, adj_p_val)

    df_integrated <- dplyr::left_join(df_integrated, da_subset, by = "Taxon")
  } else {
    df_integrated$logFC <- NA
    df_integrated$p_val <- NA
    df_integrated$adj_p_val <- NA
  }

  # Ensure order matches the plot (descending importance)
  df_integrated <- df_integrated %>%
    dplyr::arrange(match(Taxon, taxa_order))

  return(df_integrated)
}
