#' Internal Helpers for Cross-Cohort Prediction
#'
#' Small private helpers used only by \code{\link{RUMBLE_crossCohort}}. Not
#' exported.
#'
#' @name cross-cohort-internal
#' @keywords internal
#' @return No return value. This object is internal and serves only for organization.
#'
#' @importFrom stats predict
#' @importFrom dplyr bind_rows
#' @importFrom yardstick metric_set f_meas roc_auc accuracy bal_accuracy mcc
NULL

## ------------------------------------------------------------------
## Identify the outcome column inside a fold's (already CLR-transformed)
## train_data. Taxa columns are numeric (CLR log-ratios); the outcome is
## the only remaining factor column. This lets RUMBLE_crossCohort() work
## even when the source and target cohorts name their outcome variable
## differently.
## ------------------------------------------------------------------
#' @noRd
.detect_outcome_column <- function(df) {
  factor_cols <- names(df)[vapply(df, is.factor, logical(1))]
  if (length(factor_cols) != 1L) {
    stop(
      "Could not unambiguously identify the outcome column in this fold's ",
      "train_data (expected exactly one factor column after CLR ",
      "transformation; found ", length(factor_cols), "). This usually means ",
      "'fold_results' was not produced by RUMBLE(), or the object was ",
      "modified after being saved.",
      call. = FALSE
    )
  }
  factor_cols
}

## ------------------------------------------------------------------
## Align a target cohort's (raw/relative) abundance table to one source
## fold's taxa vocabulary: taxa present in the fold but absent from the
## target become an all-zero column (handled downstream by .apply_clr()'s
## pseudocount, exactly like any other structural zero); taxa present in
## the target but never seen by that fold are dropped, since the model
## has no coefficient/split for them.
## ------------------------------------------------------------------
#' @noRd
.alignTargetToVocab <- function(target_counts, vocab) {
  aligned <- matrix(
    0, nrow = nrow(target_counts), ncol = length(vocab),
    dimnames = list(rownames(target_counts), vocab)
  )
  common <- intersect(vocab, colnames(target_counts))
  if (length(common) > 0L) {
    aligned[, common] <- as.matrix(target_counts[, common, drop = FALSE])
  }
  as.data.frame(aligned)
}

## ------------------------------------------------------------------
## MCC-based weight for one (fold, model) prediction unit, mirroring the
## fallback chain already used for the algorithm-level global_weights in
## the main consensus step (MCC -> ROC AUC -> equal weight), applied here
## per (fold, model) unit instead of per algorithm -- see
## RUMBLE_crossCohort()'s documentation for why the axis changes from
## "algorithm" to "fold x algorithm".
## ------------------------------------------------------------------
#' @noRd
.mccWeightForUnit <- function(metrics_df, w_min = 0.01) {
  w <- metrics_df$.estimate[metrics_df$.metric == "mcc"]
  if (length(w) == 0L || is.na(w[1])) {
    w <- metrics_df$.estimate[metrics_df$.metric == "roc_auc"]
  }
  if (length(w) == 0L || is.na(w[1])) {
    return(1)
  }
  max(w[1], w_min)
}

## ------------------------------------------------------------------
## Combine a set of per-unit predicted probabilities (each a numeric
## vector aligned to the same target sample order) with normalized
## weights, then derive the metrics table against a common truth vector.
## Shared by the per-model and the "overall" (pooled) combination steps.
## ------------------------------------------------------------------
#' @noRd
.combineWeightedPredictions <- function(prob_list, weights, class_of_interest,
                                        negative_class, truth, sample_ids) {
  weights <- pmax(weights, 0)
  if (sum(weights) <= 0) {
    stop("All prediction units received zero weight -- cannot combine.", call. = FALSE)
  }
  weights <- weights / sum(weights)

  prob_matrix <- do.call(cbind, prob_list)
  prob_positive <- as.numeric(prob_matrix %*% weights)

  pred_class <- factor(
    ifelse(prob_positive > 0.5, as.character(class_of_interest), as.character(negative_class)),
    levels = c(as.character(negative_class), as.character(class_of_interest))
  )

  predictions <- data.frame(
    sample_id     = sample_ids,
    truth         = truth,
    prob_positive = prob_positive,
    pred_class    = pred_class,
    stringsAsFactors = FALSE
  )

  metrics_fn <- yardstick::metric_set(
    yardstick::f_meas, yardstick::accuracy,
    yardstick::bal_accuracy, yardstick::mcc
  )
  perf <- metrics_fn(predictions, truth = truth, estimate = pred_class)

  auc <- tryCatch(
    yardstick::roc_auc(predictions, truth = truth, prob_positive, event_level = "second"),
    error = function(e) NULL
  )
  if (!is.null(auc)) perf <- dplyr::bind_rows(perf, auc)

  list(predictions = predictions, metrics = perf)
}

#' Predict Across Cohorts Using Models Trained by RUMBLE()
#'
#' @description
#' Applies the models trained on a \emph{source} cohort (the per-outer-fold
#' fitted workflows saved by \code{\link{RUMBLE}} to
#' \code{<prefix>_trained_fold_models.rds}) to an entirely different,
#' \emph{target} cohort, to assess whether the biomarker signature
#' generalizes across studies (train-on-one, test-on-other).
#'
#' @details
#' For each admitted algorithm and each outer fold of the source run, the
#' fold's own fitted workflow, taxa vocabulary, and CLR pseudocount are
#' reused exactly as trained -- no re-tuning happens here. Concretely, for
#' every (fold, algorithm) combination:
#' \enumerate{
#'   \item The target cohort's raw/relative abundance table is aligned to
#'     that fold's taxa vocabulary (\code{colnames(fold$train_data)} minus
#'     the outcome column): taxa the fold expects but the target cohort
#'     does not have become an all-zero column; taxa the target cohort has
#'     but the fold never saw during training are dropped.
#'   \item CLR is applied to that aligned matrix, using \code{pseudocount}
#'     (which must match the value used to train the source models).
#'   \item \code{predict(model_fit, newdata, type = "prob")} is called --
#'     since \code{model_fit} is a fitted \pkg{workflows} object, this
#'     automatically re-applies the fold's own recipe (including
#'     \code{step_normalize()} with that fold's training mean/SD), with no
#'     separate recipe-extraction step needed.
#' }
#' The resulting per-(fold, algorithm) predicted probabilities are then
#' combined per target sample with the same MCC-based weighting scheme
#' used for the algorithm-level weighted consensus (\code{global_weights},
#' Equation 5 of the manuscript) -- but computed \emph{per (fold,
#' algorithm) unit instead of per algorithm}, using that unit's own
#' outer-fold test MCC (already stored in \code{fold_results}) as its
#' weight, floored at \code{w_min} and normalized to sum to 1. This is
#' returned both per algorithm (\code{$per_model}) -- useful to see how
#' each of RF/ENET/KNN/XGB individually transfers -- and pooled across
#' every (fold, algorithm) unit at once (\code{$overall}).
#'
#' \strong{Consistency requirements:} \code{tax_level},
#' \code{remove_unclassified}, \code{unclassified_patterns}, and
#' \code{pseudocount} must match exactly what was used to build the
#' \emph{source} cohort's \code{fold_results} -- this function has no way
#' to detect a mismatch, since that configuration is not stored inside
#' \code{fold_results} itself. A taxonomic-aggregation mismatch (e.g.
#' source aggregated to Genus, target left at ASV/OTU level) will silently
#' produce near-total vocabulary misalignment rather than an error.
#'
#' @param fold_results Either the list object saved by \code{RUMBLE()} to
#'   \code{<prefix>_trained_fold_models.rds} (one element per outer fold,
#'   each with \code{train_data} and \code{models}), or a single character
#'   path to that \code{.rds} file.
#' @param target_input A \code{phyloseq} object OR a numeric matrix/data.frame
#'   (samples as rows, taxa as columns) for the \emph{target} cohort -- same
#'   convention as \code{\link{prepareData}}'s \code{input}.
#' @param target_metadata A data.frame of target-cohort sample metadata.
#'   Required when \code{target_input} is a matrix or data.frame.
#' @param outcome_var Character. Name of the binary outcome variable in
#'   \code{target_metadata} (does not need to match the source cohort's
#'   outcome variable name).
#' @param class_of_interest Character. The level of \code{outcome_var} in
#'   the \emph{target} cohort that corresponds \emph{biologically} to the
#'   class the source models were trained to predict (e.g. both cohorts'
#'   "CRC" case group). This function does not attempt to match label
#'   strings across cohorts automatically -- \code{prob_positive} in the
#'   output is always the probability the source model assigns to its own
#'   positive class; the caller asserts the two positive classes correspond
#'   to the same biological condition.
#' @param tax_level Character or \code{NULL}. Must match the value used
#'   for the source cohort. Default \code{NULL}.
#' @param remove_unclassified,unclassified_patterns As in
#'   \code{\link{prepareData}}; must match the source cohort's settings.
#' @param pseudocount Numeric. CLR pseudocount; must match the value used
#'   to train \code{fold_results} (default \code{1e-6}, matching
#'   \code{\link{RUMBLE}}'s own default).
#' @param models Character vector of algorithm names to include (e.g.
#'   \code{c("RF", "XGB")}), or \code{NULL} (default) to use every model
#'   name found across \code{fold_results}. It is strongly recommended to
#'   pass the source run's \code{results$selected_models} here, so that
#'   only algorithms admitted by \code{metric_cutoffs} in the source cohort
#'   are used for transfer.
#' @param w_min Numeric. Minimum weight floor applied to each (fold,
#'   algorithm) unit before normalization, matching the \code{w_min = 0.01}
#'   convention already used for \code{global_weights} in the main
#'   consensus step. Default \code{0.01}.
#' @param verbose Logical. Whether to print progress messages (default
#'   \code{TRUE}).
#'
#' @return A named list:
#' \itemize{
#'   \item \code{per_model}: named list, one element per algorithm in
#'     \code{models}, each with \code{predictions} (data.frame:
#'     \code{sample_id}, \code{truth}, \code{prob_positive},
#'     \code{pred_class}), \code{metrics} (data.frame, same metric set as
#'     \code{\link{RUMBLE}}'s own fold evaluation: f_meas, accuracy,
#'     bal_accuracy, mcc, roc_auc), \code{fold_weights} (named numeric,
#'     the normalized per-fold weight actually used), and
#'     \code{taxa_alignment} (data.frame, one row per fold, reporting
#'     vocabulary size and how much of it was found in the target cohort).
#'   \item \code{overall}: same structure as one \code{per_model} entry,
#'     but pooling every (fold, algorithm) unit across all algorithms in
#'     \code{models} into a single weighted combination.
#'   \item \code{models_used}: character vector of algorithms actually
#'     used (a model requested but never found in any fold is dropped with
#'     a warning, not an error).
#'   \item \code{n_folds_total}: number of outer folds found in
#'     \code{fold_results}.
#'   \item \code{n_target_samples}: number of target-cohort samples used
#'     (after removing rows with \code{NA} in \code{outcome_var}).
#' }
#'
#' @importFrom methods is
#' @export
RUMBLE_crossCohort <- function(fold_results,
                               target_input,
                               target_metadata = NULL,
                               outcome_var,
                               class_of_interest,
                               tax_level = NULL,
                               remove_unclassified = FALSE,
                               unclassified_patterns = "uncultured|unknown|unclassified",
                               pseudocount = 1e-6,
                               models = NULL,
                               w_min = 0.01,
                               verbose = TRUE) {

  msg <- function(...) {
    if (verbose) message(...)
  }

  ## ==================================================================
  ## 1. Load / validate fold_results
  ## ==================================================================
  if (is.character(fold_results)) {
    if (length(fold_results) != 1L || !file.exists(fold_results)) {
      stop(
        "'fold_results' must be either the list object saved by RUMBLE() ",
        "to '<prefix>_trained_fold_models.rds', or a single valid path to ",
        "that .rds file.", call. = FALSE
      )
    }
    msg("Reading fold_results from: ", fold_results)
    fold_results <- readRDS(fold_results)
  }
  if (!is.list(fold_results) || length(fold_results) == 0L) {
    stop(
      "'fold_results' is empty or not a list. Expected the object saved ",
      "as '<prefix>_trained_fold_models.rds' by RUMBLE(output_dir = ...).",
      call. = FALSE
    )
  }

  n_folds_total <- length(fold_results)
  msg("Loaded ", n_folds_total, " source outer fold(s).")

  available_models <- unique(unlist(lapply(fold_results, function(fr) names(fr$models))))
  if (length(available_models) == 0L) {
    stop("No trained models found inside 'fold_results'.", call. = FALSE)
  }

  if (is.null(models)) {
    models <- available_models
    msg(
      "No 'models' specified -- using every model found across folds: ",
      paste(models, collapse = ", "),
      ". Consider passing the source run's results$selected_models to ",
      "restrict this to models admitted by metric_cutoffs."
    )
  } else {
    invalid_models <- setdiff(models, available_models)
    if (length(invalid_models) > 0L) {
      stop(
        "models not found in 'fold_results': ", paste(invalid_models, collapse = ", "),
        ". Available: ", paste(available_models, collapse = ", "), call. = FALSE
      )
    }
  }

  ## ==================================================================
  ## 2. Prepare the target cohort (same no-op min_prevalence = 0.0
  ##    pattern RUMBLE() itself uses to build analysis_df -- the actual
  ##    taxa filtering happens per-fold via vocabulary alignment below,
  ##    not via a global prevalence filter here).
  ## ==================================================================
  msg("Preparing target cohort data...")
  target_prep <- prepareData(
    input = target_input,
    metadata = target_metadata,
    tax_level = tax_level,
    min_prevalence = 0.0,
    min_abundance = 0.0001,
    remove_unclassified = remove_unclassified,
    unclassified_patterns = unclassified_patterns,
    verbose = verbose
  )
  target_counts <- target_prep$filtered_counts

  if (!outcome_var %in% colnames(target_prep$metadata)) {
    stop(
      "'", outcome_var, "' not found in target_metadata. Available: ",
      paste(colnames(target_prep$metadata), collapse = ", "), call. = FALSE
    )
  }

  target_outcome <- target_prep$metadata[[outcome_var]]
  valid_idx <- !is.na(target_outcome)
  n_before <- length(target_outcome)
  target_counts   <- target_counts[valid_idx, , drop = FALSE]
  target_outcome  <- target_outcome[valid_idx]
  if (sum(!valid_idx) > 0L) {
    msg("Removed ", sum(!valid_idx), " target samples with NA in '", outcome_var, "'")
  }

  target_outcome <- droplevels(as.factor(target_outcome))
  lvls <- levels(target_outcome)
  if (length(lvls) != 2L) {
    stop(
      "RUMBLE_crossCohort() currently supports binary outcomes only. Found ",
      length(lvls), " levels in '", outcome_var, "' (target cohort).", call. = FALSE
    )
  }
  if (!class_of_interest %in% lvls) {
    stop(
      "'", class_of_interest, "' not found in target '", outcome_var, "'. ",
      "Available classes: ", paste(lvls, collapse = ", "), call. = FALSE
    )
  }
  negative_class <- lvls[lvls != class_of_interest]
  target_outcome <- factor(target_outcome, levels = c(negative_class, class_of_interest))

  n_target_samples <- nrow(target_counts)
  target_ids <- rownames(target_counts)
  msg(
    "Target cohort ready: ", n_target_samples, " samples ('",
    negative_class, "' vs '", class_of_interest, "'), ",
    ncol(target_counts), " taxa (pre-alignment)."
  )

  ## ==================================================================
  ## 3. Per (fold, model) predictions
  ## ==================================================================
  # unit_probs / unit_weights accumulate every (fold, model) prediction
  # unit across ALL models, in the same order, so that the "overall"
  # pooled combination (section 4) can reuse them directly without
  # recomputing anything.
  unit_probs   <- list()
  unit_weights <- numeric(0)

  per_model <- list()

  for (model_name in models) {
    msg("\n--- Cross-cohort transfer: model '", model_name, "' ---")

    fold_probs   <- list()
    fold_weights <- numeric(0)
    alignment_rows <- list()

    for (fr in fold_results) {
      fold_id <- fr$fold_id
      mod_entry <- fr$models[[model_name]]
      if (is.null(mod_entry) || is.null(mod_entry$model_fit)) {
        msg("  Fold ", fold_id, ": model '", model_name, "' not available (skipped).")
        next
      }

      source_outcome_col <- .detect_outcome_column(fr$train_data)
      vocab <- setdiff(colnames(fr$train_data), source_outcome_col)

      n_matched <- length(intersect(vocab, colnames(target_counts)))
      alignment_rows[[length(alignment_rows) + 1L]] <- data.frame(
        fold = fold_id,
        n_vocab_taxa = length(vocab),
        n_matched_in_target = n_matched,
        pct_matched = round(100 * n_matched / max(length(vocab), 1L), 1)
      )
      if (n_matched / max(length(vocab), 1L) < 0.5) {
        warning(
          "Fold ", fold_id, " (model '", model_name, "'): only ", n_matched,
          " / ", length(vocab), " taxa were found in the target cohort. ",
          "Check that tax_level/remove_unclassified/unclassified_patterns ",
          "match the source cohort's settings.", call. = FALSE
        )
      }

      aligned_raw <- .alignTargetToVocab(target_counts, vocab)
      # .apply_clr() only needs *a* column matching its target_var argument
      # to know which columns are taxa vs. outcome; the values in it are
      # never used, so the target cohort's own outcome factor is fine here.
      aligned_raw[[source_outcome_col]] <- target_outcome
      aligned_clr <- .apply_clr(aligned_raw, source_outcome_col, pseudocount = pseudocount)
      newdata <- aligned_clr[, vocab, drop = FALSE]

      source_positive_level <- levels(fr$train_data[[source_outcome_col]])[2]
      prob_df  <- stats::predict(mod_entry$model_fit, newdata, type = "prob")
      prob_col <- paste0(".pred_", source_positive_level)
      if (!prob_col %in% colnames(prob_df)) {
        warning(
          "Fold ", fold_id, " (model '", model_name, "'): expected column '",
          prob_col, "' not found in predict() output (skipped).", call. = FALSE
        )
        next
      }
      prob_positive <- prob_df[[prob_col]]

      w <- .mccWeightForUnit(mod_entry$metrics, w_min = w_min)

      fold_probs[[as.character(fold_id)]] <- prob_positive
      fold_weights[as.character(fold_id)] <- w

      unit_probs[[length(unit_probs) + 1L]]   <- prob_positive
      unit_weights[length(unit_weights) + 1L] <- w
    }

    if (length(fold_probs) == 0L) {
      warning("Model '", model_name, "' was not available in any fold -- skipped entirely.", call. = FALSE)
      next
    }

    combined <- .combineWeightedPredictions(
      fold_probs, fold_weights, class_of_interest, negative_class,
      target_outcome, target_ids
    )

    per_model[[model_name]] <- list(
      predictions    = combined$predictions,
      metrics        = combined$metrics,
      fold_weights   = fold_weights / sum(fold_weights),
      taxa_alignment = dplyr::bind_rows(alignment_rows)
    )

    msg(
      "  '", model_name, "': combined ", length(fold_probs), " fold(s). ",
      "MCC = ", round(combined$metrics$.estimate[combined$metrics$.metric == "mcc"], 3),
      ", ROC AUC = ", round(combined$metrics$.estimate[combined$metrics$.metric == "roc_auc"], 3)
    )
  }

  if (length(per_model) == 0L) {
    stop("No requested model produced any usable cross-cohort prediction.", call. = FALSE)
  }

  ## ==================================================================
  ## 4. Overall (pooled across every fold x model unit)
  ## ==================================================================
  msg("\n--- Cross-cohort transfer: overall (pooled across all models/folds) ---")
  overall_combined <- .combineWeightedPredictions(
    unit_probs, unit_weights, class_of_interest, negative_class,
    target_outcome, target_ids
  )
  overall <- list(
    predictions  = overall_combined$predictions,
    metrics      = overall_combined$metrics,
    fold_weights = unit_weights / sum(unit_weights)
  )
  msg(
    "  Overall: pooled ", length(unit_probs), " (fold x model) unit(s). ",
    "MCC = ", round(overall$metrics$.estimate[overall$metrics$.metric == "mcc"], 3),
    ", ROC AUC = ", round(overall$metrics$.estimate[overall$metrics$.metric == "roc_auc"], 3)
  )

  list(
    per_model        = per_model,
    overall          = overall,
    models_used      = names(per_model),
    n_folds_total    = n_folds_total,
    n_target_samples = n_target_samples
  )
}
