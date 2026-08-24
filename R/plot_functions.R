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

  # Filter by model, if requested
  if (!is.null(model)) {
    if (!"model" %in% colnames(df)) {
      stop("Cannot filter by model: 'model' column not found in importance table.")
    }
    df <- df[df$model == model, , drop = FALSE]
  }

  # Standardize the importance column name for plotting
  if ("mean_abs_contribution" %in% colnames(df)) {
    df$importance_value <- df$mean_abs_contribution
  } else if ("contribution" %in% colnames(df)) {
    df$importance_value <- abs(df$contribution)
  } else {
    stop("Could not find a valid importance column (expected 'mean_abs_contribution').")
  }

  # Standardize the direction column (Rho)
  if ("SHAP_Rho" %in% colnames(df)) {
    df$direction <- df$SHAP_Rho
  } else if (!"direction" %in% colnames(df)) {
    df$direction <- 0
  }

  return(df)
}

#' Shorten Overly Long Taxon Labels for Display
#'
#' Internal helper that shortens full taxonomic-lineage strings -- used as
#' the \code{variable} value whenever a taxon lacks genus-level
#' classification, e.g.
#' \code{"Bacteria_Actinomycetota_Coriobacteriia_Coriobacteriales_Eggerthellaceae_"}
#' -- into a compact label suitable for axis text and facet strips.
#'
#' This is meant to be used \emph{only} via a \code{labels =} argument (e.g.
#' \code{scale_x_discrete()}/\code{scale_y_discrete()}) or a
#' \code{labeller =} argument (\code{ggplot2::as_labeller()}), never by
#' mutating the underlying \code{variable} column directly: several callers
#' (\code{plotTaxaPrevalence()}, \code{plotBiomarkerIntegrated()}) rely on
#' \code{variable} matching real column names in the abundance matrix or
#' taxon identifiers in the differential-abundance table, so the display
#' text must be able to change without touching that join key.
#'
#' @param x Character vector of taxon names.
#' @param max_width Integer. Line width (characters) used to wrap any label
#'   still too long after the lineage-collapsing step (default 24).
#' @return Character vector of the same length as \code{x}, shortened
#'   and/or wrapped for display.
#' @noRd
.shortenTaxonLabel <- function(x, max_width = 24L) {
  x <- as.character(x)

  vapply(x, function(label) {
    short <- label

    # Full unclassified lineage strings follow the
    # "Bacteria_<phylum>_..._<lowest rank>_" convention used whenever a
    # taxon has no genus-level classification. Collapse these down to just
    # the lowest (most specific) non-empty rank.
    if (grepl("^Bacteria_", label)) {
      parts <- strsplit(label, "_")[[1]]
      parts <- parts[nzchar(parts)]
      if (length(parts) >= 1) {
        short <- paste0(parts[length(parts)], " (unclassified)")
      }
    }

    if (nchar(short) > max_width) {
      short <- paste(strwrap(short, width = max_width), collapse = "\n")
    }

    short
  }, character(1), USE.NAMES = FALSE)
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

      # Compute this fold's specific ROC curve
      curve_fold <- yardstick::roc_curve(
        df_fold,
        truth = !!rlang::sym(target_var),
        !!rlang::sym(prob_col),
        event_level = "second"
      )

      # Compute this fold's specific AUC
      auc_fold <- yardstick::roc_auc(
        df_fold,
        truth = !!rlang::sym(target_var),
        !!rlang::sym(prob_col),
        event_level = "second"
      )$.estimate

      fpr_raw <- 1 - curve_fold$specificity
      tpr_raw <- curve_fold$sensitivity

      # Strict ordering to ensure the linear interpolation works correctly
      ord <- order(fpr_raw)
      fpr_raw <- fpr_raw[ord]
      tpr_raw <- tpr_raw[ord]

      # Interpolate sensitivities onto the standardized FPR grid
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

  # Aggregate results to compute the point-by-point mean and standard deviation
  df_mean_roc <- df_folds_roc %>%
    dplyr::group_by(model, fpr) %>%
    dplyr::summarize(
      mean_tpr = mean(tpr, na.rm = TRUE),
      sd_tpr = sd(tpr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Clamp the mathematical bounds so they don't exceed the biological range [0,1]
      ymin = pmax(0, mean_tpr - sd_tpr),
      ymax = pmin(1, mean_tpr + sd_tpr)
    )

  # Generate the summary metrics (Mean +/- SD) of the AUC across outer folds
  df_auc_stats <- df_auc %>%
    dplyr::group_by(model) %>%
    dplyr::summarize(
      mean_auc = mean(auc, na.rm = TRUE),
      sd_auc = sd(auc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Line break between the model name and the AUC annotation: keeps
      # each line of the facet strip short enough that it fits the strip
      # box without being clipped (strip.clip = "off" below is a second
      # safety net for any label that is still too wide).
      label = paste0(model, "\n(AUC = ", sprintf("%.3f", mean_auc), " \u00b1 ", sprintf("%.3f", sd_auc), ")")
    )

  df_mean_roc <- dplyr::left_join(df_mean_roc, df_auc_stats, by = "model")
  df_folds_roc <- dplyr::left_join(df_folds_roc, df_auc_stats, by = "model")

  p <- ggplot2::ggplot() +
    # Diagonal reference line (random classifier)
    ggplot2::geom_abline(linetype = "dashed", color = "gray60", linewidth = 0.5) +
    # Variance band (standard-deviation shading)
    ggplot2::geom_ribbon(
      data = df_mean_roc,
      ggplot2::aes(x = fpr, ymin = ymin, ymax = ymax, fill = model),
      alpha = 0.15, inherit.aes = FALSE
    ) +
    # Individual fold paths (thin, light lines in the background)
    ggplot2::geom_line(
      data = df_folds_roc,
      ggplot2::aes(x = fpr, y = tpr, group = interaction(model, fold), color = model),
      linewidth = 0.4, alpha = 0.35
    ) +
    # Consolidated mean curve for the model (thick, highlighted line)
    ggplot2::geom_line(
      data = df_mean_roc,
      ggplot2::aes(x = fpr, y = mean_tpr, color = model),
      linewidth = 1.2
    ) +
    ggplot2::facet_wrap(~ label) +
    ggsci::scale_color_nejm() +
    ggsci::scale_fill_nejm() +
    # Small margin (-0.01 to 1.01) to definitively fix visual clipping at the edges
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
      strip.text = ggplot2::element_text(face = "bold", size = 9.5, lineheight = 0.9),
      # Safety net on top of the "\n" line break in `label` above: even if
      # a label is still wider than its strip box (e.g. a long AUC SD),
      # this stops ggplot from clipping both ends of the text instead of
      # letting it overflow slightly into the panel margin.
      strip.clip = "off",
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
    # Shortens long unclassified-lineage taxon names for display only --
    # `label_with_sig` itself (used for ordering) is untouched.
    ggplot2::scale_x_discrete(labels = .shortenTaxonLabel) +
    ggplot2::scale_fill_gradient2(
      low = "#1f77b4", mid = "gray85", high = "#d62728",
      midpoint = 0, limits = c(-1, 1),
      name = "Correlation (\u03c1)"
    ) +
    ggplot2::labs(
      x = NULL,
      y = if (is.null(model)) paste0("Normalized SHAP Score (", metric_name, ")") else paste0("Raw SHAP Score (", metric_name, ")"),
      title = plot_title,
      # Line break (\n) inserted to prevent the text from being clipped
      subtitle = paste0(
        target_var, ": <- ", negative_class,
        " | ", class_of_interest, " ->\n(Direction: Spearman \u03c1 between abundance\nand SHAP impact)"
      )
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(t = 10, r = 15, b = 10, l = 10) # Safety margin
    )
}

#' Plot Confusion Matrices (Per-Fold Mean +/- SD)
#'
#' Displays confusion matrices as heatmaps for each model, one cell per
#' (Truth, Prediction) combination. Each cell shows the mean proportion of
#' that outcome computed \emph{per outer CV fold} (i.e. the confusion
#' matrix is built separately for every fold, row-normalized within that
#' fold, and then averaged), plus the standard deviation across folds and
#' the pooled total count for reference.
#'
#' This intentionally does NOT just pool every out-of-fold prediction into
#' one matrix (which was the previous behavior): a single pooled matrix
#' looks identical whether performance is stable across folds or wildly
#' inconsistent, which reads as an uninformative/"dummy" figure. Showing
#' the fold-to-fold mean +/- SD makes that variability visible, consistent
#' with how \code{plotRocCurves()} and \code{plotMetricsComparison()}
#' (with \code{fold_metrics} supplied) already report fold-level spread.
#'
#' @param final_models Named list of model results. Each element's
#'   \code{$predictions} must contain a \code{fold} column identifying the
#'   outer CV fold each row came from (as produced by the pipeline's
#'   nested-CV loop); if absent, all predictions are treated as a single
#'   fold and this degrades gracefully to the old pooled behavior (SD = 0).
#' @param target_var Character. Name of the target variable.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotConfusionMatrices <- function(final_models, target_var) {

  cm_by_fold <- purrr::imap_dfr(final_models, function(x, name) {
    preds <- x$predictions
    fold_col <- if ("fold" %in% colnames(preds)) preds$fold else rep(1L, nrow(preds))
    folds <- unique(fold_col)

    purrr::map_dfr(folds, function(f) {
      preds_fold <- preds[fold_col == f, , drop = FALSE]
      if (length(unique(preds_fold[[target_var]])) < 2 || nrow(preds_fold) < 2) {
        return(NULL)
      }
      cm <- yardstick::conf_mat(
        preds_fold,
        truth = !!rlang::sym(target_var),
        estimate = .pred_class
      )
      df <- as.data.frame(cm$table)
      df$model <- name
      df$fold <- f
      df
    })
  })

  if (nrow(cm_by_fold) == 0) {
    return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No confusion matrix data available"))
  }

  # Row-normalize WITHIN each fold (proportion of each true class predicted
  # as each class, for that fold) -- the fold-level counterpart of the
  # pooled `Freq / sum(Freq) * 100` used previously.
  cm_by_fold <- cm_by_fold %>%
    dplyr::group_by(model, fold, Truth) %>%
    dplyr::mutate(fold_total = sum(Freq)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Prop = ifelse(fold_total > 0, Freq / fold_total * 100, NA_real_))

  cm_summary <- cm_by_fold %>%
    dplyr::group_by(model, Truth, Prediction) %>%
    dplyr::summarize(
      mean_prop  = mean(Prop, na.rm = TRUE),
      sd_prop    = stats::sd(Prop, na.rm = TRUE),
      total_freq = sum(Freq),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      sd_prop  = ifelse(is.na(sd_prop), 0, sd_prop),
      text_col = ifelse(mean_prop > 50, "white", "black")
    )

  ggplot2::ggplot(
    cm_summary,
    ggplot2::aes(x = Prediction, y = Truth, fill = mean_prop)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = paste0(
          sprintf("%.1f", mean_prop), "% ± ", sprintf("%.1f", sd_prop), "%\n",
          "(", total_freq, " total)"
        ),
        color = text_col
      ),
      size = 3.3, fontface = "bold", lineheight = 0.95
    ) +
    ggplot2::facet_wrap(~ model) +
    ggplot2::scale_fill_gradient(
      low = "#eff3ff",
      high = "#08519c",
      limits = c(0, 100),
      name = "Mean Proportion (%)"
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      title = "Out-of-Fold Confusion Matrices (Per-Fold Mean ± SD)",
      subtitle = "Each cell: mean proportion across outer CV folds ± SD (pooled total count in parentheses)",
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
#' Produces a grouped barplot comparing performance metrics across all
#' models. When \code{fold_metrics} is supplied, bars show the mean across
#' outer CV folds with an error bar for \eqn{\pm}1 SD -- consistent with the
#' fold-level variability already shown in \code{plotRocCurves()} -- instead
#' of a single pooled-OOF point estimate with no indication of spread.
#'
#' @param final_models Named list of model results. Used only as the
#'   fallback data source when \code{fold_metrics} is not supplied (kept for
#'   backward compatibility with existing calls).
#' @param fold_metrics Optional data.frame with one row per
#'   model/fold/metric (columns \code{model}, \code{fold}, \code{.metric},
#'   \code{.estimate}), such as the \code{all_metrics} table already
#'   assembled internally by \code{RUMBLE()}'s outer-CV loop. When supplied
#'   (recommended), bars show mean \eqn{\pm} SD across folds. When
#'   \code{NULL} (default), falls back to the single pooled-OOF point
#'   estimate in \code{final_models[[.]]$metrics}, with no error bars.
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotMetricsComparison <- function(final_models, fold_metrics = NULL) {

  has_folds <- !is.null(fold_metrics) && nrow(fold_metrics) > 0

  if (has_folds) {
    plot_df <- fold_metrics %>%
      dplyr::group_by(model, .metric) %>%
      dplyr::summarize(
        Value = mean(.estimate, na.rm = TRUE),
        sd_estimate = stats::sd(.estimate, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::rename(Metric = .metric) %>%
      dplyr::mutate(
        sd_estimate = ifelse(is.na(sd_estimate), 0, sd_estimate),
        ymin = Value - sd_estimate,
        ymax = Value + sd_estimate
      )
  } else {
    metrics_df <- purrr::imap_dfr(final_models, function(obj, name) {
      obj$metrics$model <- name
      obj$metrics
    })

    plot_df <- metrics_df %>%
      dplyr::select(model, .metric, .estimate) %>%
      tidyr::pivot_wider(
        names_from = .metric, values_from = .estimate
      ) %>%
      tidyr::pivot_longer(
        cols = -model, names_to = "Metric",
        values_to = "Value"
      )
  }

  # Metrics like `mcc` range [-1, 1] and can legitimately go negative (e.g.
  # a near-random model), unlike the [0, 1] metrics plotted alongside it.
  # A fixed `ylim(0, 1)` would silently DROP those bars instead of just
  # constraining the visible area, so the lower bound is only relaxed
  # below 0 when the data actually requires it.
  lower_bound <- min(0, min(if (has_folds) plot_df$ymin else plot_df$Value, na.rm = TRUE))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = model, y = Value, fill = Metric)
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.7
    )

  if (has_folds) {
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ymin, ymax = ymax),
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.2, linewidth = 0.4, color = "gray30"
    )

    # Individual outer-fold estimates, one dot per (model, metric, fold),
    # drawn on top of the mean bar +/- SD error bar above. A bar with an
    # error bar alone can still read as a synthetic/"dummy" summary; seeing
    # the actual per-fold values scattered underneath (the same idea
    # already used for the thin per-fold lines in plotRocCurves()) is what
    # makes it visibly real fold-level data. Positions are computed
    # manually (not via position_jitterdodge()) so placement is fully
    # deterministic -- each fold gets a fixed horizontal slot based on its
    # rank within the (model, metric) group, fanned out symmetrically
    # around the bar's dodge position, rather than relying on RNG-based
    # jitter (which would also require save/restore of the caller's RNG
    # state, as done in plotShapBeeswarm()).
    model_levels  <- sort(unique(as.character(plot_df$model)))
    metric_levels <- sort(unique(as.character(plot_df$Metric)))
    n_metric <- length(metric_levels)
    dodge_width <- 0.8
    slot_width  <- dodge_width / n_metric

    fold_points <- fold_metrics %>%
      dplyr::rename(Metric = .metric, Value = .estimate) %>%
      dplyr::group_by(model, Metric) %>%
      dplyr::mutate(
        .n_folds    = dplyr::n_distinct(fold),
        .fold_rank  = dplyr::dense_rank(fold) - (.n_folds + 1) / 2
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        .model_idx  = match(model, model_levels),
        .metric_idx = match(Metric, metric_levels),
        x_point     = .model_idx +
          (.metric_idx - (n_metric + 1) / 2) * slot_width +
          .fold_rank * (slot_width / (pmax(.n_folds, 1) + 1))
      )

    p <- p + ggplot2::geom_point(
      data = fold_points,
      ggplot2::aes(x = x_point, y = Value),
      inherit.aes = FALSE,
      size = 1.3, alpha = 0.6, shape = 16, color = "black"
    )
  }

  p +
    ggsci::scale_fill_nejm() +
    ggplot2::labs(
      title = "Model Performance Comparison",
      subtitle = if (has_folds) {
        # Kept short on purpose: ggplot2 does not auto-wrap plot subtitles,
        # so a longer version of this line was overflowing past the right
        # edge of the saved PNG at the figure's normal width.
        "Bars: fold mean | Error bars: ± 1 SD | Dots: each outer fold"
      } else {
        NULL
      },
      x = "Model", y = "Metric Value", fill = "Metric"
    ) +
    ggplot2::theme_classic(base_size = 14) +
    # coord_cartesian() only zooms the visible area -- unlike ylim()/
    # scale_y_continuous(limits=...), it never drops data outside range.
    ggplot2::coord_cartesian(ylim = c(lower_bound, 1)) +
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
plotShapBeeswarm <- function(shap_raw, prediction_data = NULL, target_var = NULL, top_n = 20L,
                             model = NULL) {

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

  # Get the top features ranked by mean absolute SHAP
  top_features <- shap_raw %>%
    dplyr::group_by(variable) %>%
    dplyr::summarize(mean_abs = mean(abs(contribution), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abs)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(variable)

  shap_filtered <- shap_raw %>%
    dplyr::filter(variable %in% top_features)

  # Normalize feature values (Min-Max scaling per taxon)
  shap_filtered <- shap_filtered %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(
      feature_val_scaled = (feature_value - min(feature_value, na.rm = TRUE)) /
        (max(feature_value, na.rm = TRUE) - min(feature_value, na.rm = TRUE) + 1e-6)
    ) %>%
    dplyr::ungroup()

  # Join with the true classes
  has_classes <- !is.null(prediction_data) && !is.null(target_var)
  if (has_classes) {
    class_df <- data.frame(
      observation_id = seq_len(nrow(prediction_data)),
      True_Class = prediction_data[[target_var]]
    )
    shap_filtered <- dplyr::left_join(shap_filtered, class_df, by = "observation_id")
  }

  # NOTE: an earlier version of this function subsampled very dense
  # features (e.g. many SHAP repetitions per observation) down to a fixed
  # point cap before plotting, to keep the swarm legible. That was
  # reverted: for a SHAP beeswarm, every observation's value is part of
  # the actual result, and silently dropping a random subset of them
  # (even for a repeated-estimation dataset where the dropped points are
  # near-duplicates of ones that are kept) misrepresents the underlying
  # data. Visual crowding on very dense features is handled below purely
  # through rendering (density-scaled transparency + `corral = "random"`),
  # which never removes a single point.
  #
  # `.point_alpha` scales down per-point opacity as the number of points
  # sharing a (variable, model) row grows, so that on a sparse feature
  # points stay clearly visible, while on a feature with thousands of
  # near-duplicate repeated-SHAP points the overlap blends into a smooth
  # density-like gradient instead of one flat opaque block -- every point
  # is still drawn, just more transparently. Capped at 1 (not just <1) so
  # ordinary feature sizes (a few hundred points or fewer -- the common
  # case with repetitions = 1) stay fully opaque, unchanged from before;
  # the fade-out only kicks in once a feature's point count climbs well
  # past that.
  shap_filtered <- shap_filtered %>%
    dplyr::group_by(variable, model) %>%
    dplyr::mutate(
      .point_alpha = pmax(0.05, pmin(1, 300 / dplyr::n()))
    ) %>%
    dplyr::ungroup()

  plot_title <- if (is.null(model)) {
    "SHAP Distribution by Model"
  } else {
    paste0("SHAP Distribution: ", model)
  }

  # Base plot
  p <- ggplot2::ggplot(
    shap_filtered,
    ggplot2::aes(
      x = contribution,
      y = reorder(variable, abs(contribution), FUN = mean)
    )
  ) +
    # Central zero line
    ggplot2::geom_vline(xintercept = 0, color = "black")

  # Add the points using the swarm algorithm. `corral = "random"` (rather
  # than "gutter") jitters points that don't fit the swarm algorithm's
  # exact placement randomly within the row instead of piling them up
  # solid against the row boundary -- "gutter" is what produced the dense
  # rectangular mass covering the whole plot area on datasets with many
  # repeated points per feature. `alpha = .point_alpha` (computed above,
  # scaled by how many points share that row) is what keeps a very dense
  # feature from reading as one flat opaque block, WITHOUT dropping any
  # point: every row of `shap_filtered` is still passed to the geom.
  if (has_classes) {
    p <- p + ggbeeswarm::geom_beeswarm(
      ggplot2::aes(color = feature_val_scaled, shape = True_Class, alpha = .point_alpha),
      size = 1,
      method = "swarm",
      corral = "random",
      corral.width = 0.85,
      priority = "density",
      cex = 0.35
    ) +
      ggplot2::labs(shape = "Class")
  } else {
    p <- p + ggbeeswarm::geom_beeswarm(
      ggplot2::aes(color = feature_val_scaled, alpha = .point_alpha),
      size = 1,
      method = "swarm",
      corral = "random",
      corral.width = 0.85,
      priority = "density",
      cex = 0.35
    )
  }

  # Apply the Viridis gradient and classic theme. `scale_alpha_identity()`
  # tells ggplot to use `.point_alpha`'s numeric values directly as opacity
  # (no separate alpha legend) -- the color scale no longer bakes in its
  # own fixed alpha, since per-point density-scaled alpha now does that job.
  p <- p + viridis::scale_color_viridis(
    option = "H",
    breaks = c(0, 1),
    labels = c("Low", "High"),
    name = "Feature Value"
  ) +
    ggplot2::scale_alpha_identity() +
    ggplot2::labs(
      x = "SHAP Value (Impact on Model Output)",
      y = "Feature",
      title = plot_title,
      subtitle = "Point opacity scales down with the number of (near-)duplicate values per feature"
    ) +
    # Shortens long unclassified-lineage taxon names for display only --
    # the underlying `variable` values (used for ordering/grouping above)
    # are untouched.
    ggplot2::scale_y_discrete(labels = .shortenTaxonLabel) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "italic")
    )

  # If this is the consensus plot, split it into panels
  if (is.null(model) && length(unique(shap_filtered$model)) > 1) {
    p <- p + ggplot2::facet_wrap(~ model, scales = "free_x")
  }

  return(p)
}

#' Plot SHAP Dependence
#'
#' Displays SHAP dependence plots for the top features: feature abundance
#' (CLR-transformed) on the x-axis against the raw SHAP contribution on the
#' y-axis, one panel per taxon. Complements \code{plotShapBeeswarm()}: the
#' beeswarm only insinuates the shape of the feature-to-SHAP relationship
#' through a color gradient, while this plot shows it directly as a
#' scatter, making threshold effects and non-monotone patterns easier to
#' inspect -- relevant because Spearman's rho (Section 2.4.3 of the
#' manuscript) is uninformative precisely for non-monotone relationships.
#'
#' @param shap_raw A data.frame of raw SHAP values (must contain
#'   \code{variable}, \code{contribution}, \code{feature_value}, and
#'   \code{model} columns, as returned by \code{computeFeatureImportance()}).
#' @param top_n Integer. Number of top features to show (default 20).
#' @param model Character or \code{NULL}. Optional model name to isolate a
#'   specific model in the dependence plot.
#' @param ncol Integer. Number of facet columns (default 4).
#'
#' @return A \code{ggplot} object.
#'
#' @export
plotShapDependence <- function(shap_raw, top_n = 20L, model = NULL, ncol = 4L) {

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

  # Get the top features ranked by mean absolute SHAP (same ranking logic
  # as plotShapBeeswarm(), so the two plots show the same feature set).
  top_features <- shap_raw %>%
    dplyr::group_by(variable) %>%
    dplyr::summarize(mean_abs = mean(abs(contribution), na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_abs)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(variable)

  shap_filtered <- shap_raw %>%
    dplyr::filter(variable %in% top_features)

  plot_title <- if (is.null(model)) {
    "SHAP Dependence (Consensus Models)"
  } else {
    paste0("SHAP Dependence: ", model)
  }

  p <- ggplot2::ggplot(
    shap_filtered,
    ggplot2::aes(x = feature_value, y = contribution)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    ggplot2::geom_point(alpha = 0.5, size = 1, color = "steelblue") +
    # `labeller` shortens long unclassified-lineage taxon names for
    # display only -- the `variable` column driving the facet split is
    # untouched. Without this, those names overflow the facet strip box
    # and get clipped on both ends (ggplot2 centers strip text).
    ggplot2::facet_wrap(~ variable, scales = "free", ncol = ncol,
                        labeller = ggplot2::as_labeller(.shortenTaxonLabel)) +
    ggplot2::labs(
      x = "Feature Value (CLR-transformed abundance)",
      y = "SHAP Value (Impact on Model Output)",
      title = plot_title
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "italic", size = 8.5, lineheight = 0.9),
      # Safety net: even a shortened label can still overflow a very
      # narrow facet box (many columns); this avoids clipping either end.
      strip.clip = "off"
    )

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
    # Shortens long unclassified-lineage taxon names for display only --
    # `variable` (used for ordering above) is untouched.
    ggplot2::scale_y_discrete(labels = .shortenTaxonLabel) +
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
  # 1. Prepare data
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
  ## Panel A: ML Importance
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
    # Shortens long unclassified-lineage taxon names for display only --
    # `label_with_sig` (used for ordering/joins across panels B-D) is
    # untouched.
    ggplot2::scale_y_discrete(labels = .shortenTaxonLabel) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9, face = "italic"),
                   panel.grid.major.y = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9, color = "gray30"))

  ## ==================================================================
  ## Panel B: Taxa Prevalence
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
  ## Panels C & D: Differential Abundance (Reintegrated)
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
#' test folds, separated by true class. This reveals class separability
#' (discriminative power) that aggregate metrics cannot capture; it does
#' not assess probability calibration (the agreement between predicted
#' probabilities and observed frequencies), which this function does not
#' compute.
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

  # Aggregate the permutation values (averaging across models)
  perm_agg <- permutation_top %>%
    dplyr::group_by(variable) %>%
    dplyr::summarize(
      perm_importance = mean(mean_delta_loss, na.rm = TRUE),
      .groups = "drop"
    )

  # Get the global SHAP
  shap_df <- global_importance %>%
    dplyr::select(variable, mean_abs_contribution) %>%
    dplyr::rename(shap_importance = mean_abs_contribution)

  # Join the two metrics
  df_merged <- dplyr::inner_join(shap_df, perm_agg, by = "variable")

  if (nrow(df_merged) == 0) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::labs(title = "No overlapping features for SHAP vs Permutation"))

  # Define the medians used to draw the quadrants
  med_shap <- stats::median(df_merged$shap_importance, na.rm = TRUE)
  med_perm <- stats::median(df_merged$perm_importance, na.rm = TRUE)

  # Flag whoever is in the upper-right quadrant (the most robust ones)
  df_merged <- df_merged %>%
    dplyr::mutate(
      is_top_quadrant = (shap_importance > med_shap) & (perm_importance > med_perm),
      point_color = ifelse(is_top_quadrant, "#d62728", "gray60"),
      # Shortened for display only; `variable` itself is untouched.
      label = ifelse(is_top_quadrant, .shortenTaxonLabel(variable), "") # Only label the top ones
    )

  p <- ggplot2::ggplot(df_merged, ggplot2::aes(x = perm_importance, y = shap_importance)) +
    # Quadrant lines
    ggplot2::geom_vline(xintercept = med_perm, linetype = "dashed", color = "gray80") +
    ggplot2::geom_hline(yintercept = med_shap, linetype = "dashed", color = "gray80") +
    # Fixed points (highlighting the top ones in red)
    ggplot2::geom_point(ggplot2::aes(color = point_color), size = 3, alpha = 0.8) +
    ggplot2::scale_color_identity() +
    # Trend line
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "gray30", linetype = "dashed", linewidth = 0.5) +
    ggplot2::labs(
      title = "SHAP vs Permutation Importance",
      subtitle = "Highlighting robust biomarkers (Upper-Right Quadrant)",
      x = "Permutation Importance (Mean Delta Loss)",
      y = "SHAP Importance (Mean |SHAP|)"
    ) +
    ggplot2::theme_minimal(base_size = 12)

  # Use ggrepel if installed (avoids overlapping text labels).
  # `max.overlaps = Inf` forces every label to be drawn even when there is
  # no room to place it without overlapping another -- with few labeled
  # points (small/test datasets) this is harmless, but with a larger
  # consensus feature set (real data) it is precisely what causes taxon
  # names to pile up on top of each other. Capping it lets ggrepel drop
  # (rather than force-overlap) whichever labels still don't fit after
  # `max.time`/`max.iter` of repel attempts.
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      ggplot2::aes(label = label),
      size = 3.5,
      fontface = "italic",
      box.padding = 0.5,
      max.overlaps = 15,
      min.segment.length = 0,
      max.time = 2,
      max.iter = 20000,
      seed = 42,
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
