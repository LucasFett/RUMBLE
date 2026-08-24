#' Create DALEX Explainer Objects
#'
#' Wraps each fitted tidymodels workflow into a DALEX explainer object.
#'
#' @param final_models List of model results from the pipeline. Each element
#'   must contain a \code{model_fit} component (the fitted workflow).
#' @param data Training data.frame (features + target).
#' @param target_var Name of the target variable.
#' @param class_of_interest Character. The class of interest for SHAP values.
#' @param verbose Logical. Whether to print progress messages (default TRUE).
#'
#' @return A list of DALEX explainer objects, one for each model.
#'
#' @importFrom DALEX explain
#' @importFrom purrr imap
#' @importFrom stats predict
#' @export
createExplainers <- function(final_models, data, target_var,
                             class_of_interest, verbose = TRUE) {
  if (verbose) message("Creating DALEX explainers...")

  # Note: If parallelism is enabled at the importance-computation step
  # (via furrr/future), XGBoost models can trigger a serialization error
  # (an externalptr cannot be passed to multisession workers without
  # explicit serialization). If this happens, consider using
  # future::plan(sequential) for this step only, or serializing the XGB
  # model via xgb.save.raw()/xgb.load.raw().

  # Validate class_of_interest
  if (missing(class_of_interest)) {
    stop("'class_of_interest' is required in createExplainers().")
  }

  # Ensure class_of_interest is the second level (positive class for DALEX)
  lvls <- levels(data[[target_var]])
  if (lvls[2L] != class_of_interest) {
    stop("Expected '", class_of_interest, "' to be the second factor level, ",
         "but found '", lvls[2L], "'. Please ensure factor levels are reordered ",
         "before calling createExplainers().")
  }

  positive_class <- class_of_interest

  # Convert target to numeric 0/1 for DALEX
  y_numeric <- ifelse(
    data[[target_var]] == positive_class, 1L, 0L
  )

  # Remove target from features. (Note: the previous version also
  # excluded a "case_weight" column here, a remnant of the
  # case-weights-via-data-column approach. That approach was replaced by
  # native per-engine class_weights (see .buildWorkflows()), which does
  # not attach any extra column to `data` -- so the exclusion was
  # simplified back to just target_var.)
  X <- data[, !colnames(data) %in% c(target_var), drop = FALSE]

  purrr::imap(final_models, function(model_obj, model_name) {
    fitted_wf <- model_obj$model_fit

    DALEX::explain(
      model = fitted_wf,
      data = X,
      y = y_numeric,
      label = model_name,
      type = "classification",
      verbose = FALSE,
      predict_function = function(m, newdata) {
        pred <- stats::predict(m, newdata, type = "prob")
        pred[[paste0(".pred_", positive_class)]]
      }
    )
  })
}


#' Compute Feature Importance via SHAP and Permutation
#'
#' Calculates SHAP values and permutation importance using parallel processing.
#' Supports both exact SHAP (via DALEX) and fast approximation (via fastshap).
#'
#' \strong{Methodological Note:} SHAP values are returned as raw (unscaled)
#' contributions. No per-fold normalization is applied to avoid distorting the
#' consensus. Direction is computed as a descriptive metric (mean SHAP for
#' observations above the feature median), NOT as a hypothesis test.
#'
#' @param explainers List of DALEX explainers (output of \code{createExplainers}).
#' @param data Data.frame for SHAP predictions (can be train or test).
#' @param target_var Name of the target variable.
#' @param class_of_interest Character. The class of interest.
#' @param top_n Number of top features to retain (default 20).
#' @param repetitions Number of repetitions (B) for importance calculation
#'   (default 10).
#' @param n_cores Number of cores for parallel processing (default 1).
#' @param shap_method Character. Method for SHAP calculation: \code{"exact"}
#'   or \code{"fast"}.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{permutation_raw}: Raw permutation importance data.
#'   \item \code{permutation_top}: Top features by permutation importance.
#'   \item \code{shap_raw}: Raw SHAP values (unscaled).
#'   \item \code{shap_top}: Top features by SHAP values per model.
#' }
#'
#' @importFrom DALEX model_parts predict_parts
#' @importFrom dplyr group_by summarize arrange slice_head ungroup desc
#'      n_distinct bind_rows mutate
#' @importFrom furrr future_map_dfr furrr_options
#' @importFrom future plan multisession sequential
#' @importFrom purrr map_dfr imap
#' @importFrom tidyr pivot_longer
#' @importFrom stats cor median
#' @export
computeFeatureImportance <- function(explainers,
                                     data,
                                     target_var,
                                     class_of_interest,
                                     top_n = 20L,
                                     repetitions = 10L,
                                     n_cores = 1L,
                                     shap_method = "exact",
                                     verbose = TRUE) {

  msg <- function(...) {
    if (verbose) message(...)
  }

  ## Validate shap_method
  shap_method <- match.arg(shap_method, c("exact", "fast"))

  msg("------------------------------------------------------")
  msg("Feature importance analysis (RUMBLE)")
  msg("Models: ", length(explainers),
      " | Top N: ", top_n,
      " | Repetitions: ", repetitions,
      " | SHAP method: ", shap_method)
  msg("------------------------------------------------------")

  # --- PARALLEL SETUP ---
  # The XGBoost booster ('xgb.Booster') holds its state as a C pointer
  # ('externalptr'), which CANNOT be serialized to future/furrr
  # (multisession) worker processes. This has already caused a real crash
  # during manual testing: the worker that received the XGB model had its
  # connection broken ('error writing to connection') and became an
  # orphan process, because once the connection breaks, not even
  # future::plan("sequential") can shut it down cleanly -- the zombie
  # process survives and can interfere with parallel clusters opened by
  # later RUMBLE() calls (via .tuneModels()/parallel::makeCluster()).
  #
  # Rather than disabling multisession for the ENTIRE model set whenever
  # an XGBoost explainer is present (the previous, coarser strategy),
  # explainers are split into an "xgb" group (always processed
  # sequentially, never touching the future/furrr backend) and an
  # "other" group (RF/ENET/KNN, eligible for multisession when
  # n_cores > 1). This keeps the preventive guarantee for XGBoost while
  # no longer forcing RF/ENET/KNN into single-core execution just
  # because XGBoost happens to be in the same model pool.
  is_xgb <- vapply(explainers, function(e) identical(e$label, "XGB"), logical(1))
  xgb_explainers   <- explainers[is_xgb]
  other_explainers <- explainers[!is_xgb]

  effective_n_cores <- if (length(other_explainers) > 0) n_cores else 1L

  if (effective_n_cores > 1) {
    msg("Initializing parallel backend with ", effective_n_cores, " cores for ",
        length(other_explainers), " non-XGBoost model(s)...")
    future::plan(future::multisession, workers = effective_n_cores)
  } else {
    future::plan(future::sequential)
  }
  if (length(xgb_explainers) > 0) {
    msg(length(xgb_explainers), " XGBoost model(s) detected: the booster's ",
        "externalptr cannot be serialized to multisession workers, so ",
        "XGBoost is always processed sequentially, isolated from the ",
        "parallel backend above (RF/ENET/KNN above remain eligible for ",
        "multisession when n_cores > 1) -- tuning/fitting parallelism ",
        "elsewhere in the pipeline is unaffected.")
  }
  # on.exit() ensures the plan reverts to sequential even if the
  # SHAP/permutation importance computation below throws an error. This
  # complements the on.exit() already present at the top of RUMBLE():
  # that one resets the session's plan at the end of the call, but it
  # does not cleanly shut down a worker whose connection has already
  # broken -- which is why routing XGBoost around the multisession
  # backend entirely (above) is the primary defense, and this is a
  # second layer of safety.
  on.exit(future::plan(future::sequential), add = TRUE)

  # --- MAPPER DEFINITIONS ---
  # run_map(): used only for the "other" (non-XGBoost) group. Defensive
  # fallback: if parallel execution still fails (for example, due to
  # some other non-exportable object besides the XGBoost case already
  # routed around above), re-run the step sequentially with a warning,
  # instead of aborting the entire pipeline.
  run_map <- function(x, fn, seed = TRUE) {
    if (effective_n_cores > 1) {
      tryCatch({
        furrr::future_map_dfr(x, fn,
                              .options = furrr::furrr_options(seed = seed))
      }, error = function(e) {
        warning(
          "Parallel execution failed (", conditionMessage(e), "). ",
          "Falling back to sequential execution for this step.",
          call. = FALSE
        )
        purrr::map_dfr(x, fn)
      })
    } else {
      purrr::map_dfr(x, fn)
    }
  }
  # run_map_seq(): used only for the XGBoost group -- always plain
  # purrr::map_dfr(), never touches future/furrr, regardless of
  # effective_n_cores (which may be > 1 because of the "other" group).
  run_map_seq <- function(x, fn, seed = TRUE) purrr::map_dfr(x, fn)

  # --- PREPARE DATA ---
  X_data <- data[, !colnames(data) %in% c(target_var), drop = FALSE]
  y_numeric <- ifelse(data[[target_var]] == class_of_interest, 1L, 0L)

  # --- 1. PERMUTATION IMPORTANCE ---
  msg("Computing permutation importance...")

  compute_perm_for_model <- function(exp) {
    parts <- DALEX::model_parts(
      exp,
      data = X_data,
      y = y_numeric,
      type = "variable_importance",
      B = repetitions
    )
    parts$model <- exp$label
    parts
  }

  perm_other <- if (length(other_explainers) > 0) {
    run_map(other_explainers, compute_perm_for_model)
  } else {
    data.frame()
  }
  perm_xgb <- if (length(xgb_explainers) > 0) {
    run_map_seq(xgb_explainers, compute_perm_for_model)
  } else {
    data.frame()
  }
  perm <- dplyr::bind_rows(perm_other, perm_xgb)



  # --- 2. SHAP VALUES ---
  if (shap_method == "exact") {
    msg("Computing SHAP values for all observations (EXACT method)...")

    compute_shap_for_model <- function(exp, mapper) {
      msg(paste0("  -> Extracting SHAP for model: ", exp$label))

      mapper(seq_len(nrow(X_data)), function(i) {
        single_obs <- X_data[i, , drop = FALSE]

        shap_val <- DALEX::predict_parts(
          exp,
          new_observation = single_obs,
          type = "shap",
          B = repetitions
        )

        shap_val$model <- exp$label
        shap_val$feature_value <- as.numeric(sub(".*= ", "", as.character(shap_val$variable)))
        shap_val$observation_id <- i

        return(shap_val)
      })
    }

    # Same split as the permutation step above: the per-observation loop
    # for non-XGBoost models may use the multisession backend (run_map);
    # XGBoost's per-observation loop always uses run_map_seq(), never
    # touching future/furrr, regardless of effective_n_cores.
    shap_other <- if (length(other_explainers) > 0) {
      purrr::map_dfr(other_explainers, function(exp) compute_shap_for_model(exp, run_map))
    } else {
      data.frame()
    }
    shap_xgb <- if (length(xgb_explainers) > 0) {
      purrr::map_dfr(xgb_explainers, function(exp) compute_shap_for_model(exp, run_map_seq))
    } else {
      data.frame()
    }
    shap <- dplyr::bind_rows(shap_other, shap_xgb)

    # Clean variable names
    shap$variable <- sub(" =.*", "", as.character(shap$variable))

  } else if (shap_method == "fast") {
    msg("Computing SHAP values using FAST approximation method...")

    if (!requireNamespace("fastshap", quietly = TRUE)) {
      stop("The 'fastshap' package is required for shap_method='fast'. ",
           "Install it with: install.packages('fastshap')")
    }

    shap <- purrr::map_dfr(explainers, function(exp) {
      msg(paste0("  -> Computing fast SHAP for model: ", exp$label))

      fitted_model <- exp$model
      X_pred <- X_data

      pred_fn <- function(object, newdata) {
        stats::predict(object, newdata, type = "prob")[[paste0(".pred_", class_of_interest)]]
      }

      shap_vals <- fastshap::explain(
        object = fitted_model,
        X = exp$data,         # The background (training set saved by DALEX)
        newdata = X_pred,     # The current sample being explained (test or train)
        pred_wrapper = pred_fn,
        nsim = repetitions
      )

      shap_long <- tidyr::pivot_longer(
        data.frame(observation_id = seq_len(nrow(shap_vals)), shap_vals),
        cols = -observation_id,
        names_to = "variable",
        values_to = "contribution"
      )

      X_long <- X_pred %>%
        dplyr::mutate(observation_id = dplyr::row_number()) %>%
        tidyr::pivot_longer(
          cols = -observation_id,
          names_to = "variable",
          values_to = "feature_value"
        )

      shap_long <- shap_long %>%
        dplyr::left_join(X_long, by = c("observation_id", "variable")) %>%
        dplyr::mutate(model = exp$label)

      return(shap_long)
    })
  }

  # Parallel backend cleanup is now handled via on.exit() in the PARALLEL
  # SETUP block above -- it covers both the success and error paths.

  # --- 3. AGGREGATE RESULTS ---
  msg("Aggregating results (raw SHAP values, no per-fold normalization)...")

  ## Top permutation
  top_perm <- perm %>%
    dplyr::group_by(model) %>%
    dplyr::mutate(
      full_model_loss = mean(dropout_loss[variable == "_full_model_"], na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!variable %in% c("_full_model_", "_baseline_")) %>%
    dplyr::mutate(delta_loss = dropout_loss - full_model_loss) %>%
    dplyr::group_by(model, variable) %>%
    dplyr::summarize(mean_delta_loss = mean(delta_loss, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(model, dplyr::desc(mean_delta_loss)) %>%
    dplyr::group_by(model) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()

  ## --- 3.1 Direction Metric (Descriptive, NOT a hypothesis test) ---
  ## For each model-variable pair, compute direction as the mean SHAP

  ## contribution for observations where feature_value > median.
  ## This is biologically interpretable: positive direction means
  ## "higher abundance -> higher probability of class_of_interest"

  ## Top SHAP per model (using RAW contributions)
  top_shap <- shap %>%
    dplyr::group_by(model, variable) %>%
    dplyr::summarize(
      mean_abs_contribution = mean(abs(contribution), na.rm = TRUE),
      sd_abs_contribution = sd(abs(contribution), na.rm = TRUE),
      SHAP_Rho = .compute_direction_internal(feature_value, contribution),
      .groups = "drop"
    ) %>%
    dplyr::group_by(model) %>%
    dplyr::arrange(model, dplyr::desc(mean_abs_contribution)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()


  msg("Feature importance analysis complete.")

  list(
    permutation_raw = perm,
    permutation_top = top_perm,
    shap_raw = shap,
    shap_top = list(spearman = top_shap)
  )
}

#' Internal Direction Computation (Spearman Correlation)
#'
#' Computes the direction metric as the Spearman correlation between feature
#' values and SHAP contributions. This captures the monotonic relationship
#' between feature abundance and model prediction impact, respecting the
#' biological transition from absence (0) to hyperproliferation.
#'
#' @param feature_value Numeric vector of feature values (abundances).
#' @param contribution Numeric vector of SHAP contributions.
#' @param weight Optional numeric vector of per-observation weights (e.g.,
#'   the per-model MCC-derived consensus weight from Equation 5, recycled
#'   across that model's rows). When \code{NULL} (default), an unweighted
#'   Spearman correlation is returned -- unchanged behavior, correct for
#'   per-model aggregates where weight is constant within the group. When
#'   supplied, a weighted Spearman is computed as a weighted Pearson
#'   correlation on the (unweighted) ranks of \code{feature_value} and
#'   \code{contribution}, so that observations coming from higher-weight
#'   (higher-MCC) models contribute proportionally more to the consensus
#'   direction estimate -- consistent with \code{mean_abs_contribution},
#'   which already uses \code{stats::weighted.mean(..., w = weight)} for
#'   the magnitude metric in the same consensus step.
#' @return A single numeric value (Spearman rho, ranging from -1 to 1).
#'   Returns 0 if correlation cannot be computed (e.g., zero variance).
#' @noRd
.compute_direction_internal <- function(feature_value, contribution, weight = NULL) {
  # If the variable is constant (e.g., 100% zeros in the fold), the variance is 0
  if (length(unique(feature_value)) < 2 || length(unique(contribution)) < 2) {
    return(0)
  }

  if (is.null(weight)) {
    # Spearman captures non-linearity and handles the CLR's zero-tail well
    rho <- suppressWarnings(stats::cor(feature_value, contribution, method = "spearman"))
    if (is.na(rho)) return(0)
    return(rho)
  }

  # Weighted Spearman: (unweighted) ranks submitted to a weighted Pearson
  # correlation -- standard convention for weighted Spearman.
  rx <- rank(feature_value, ties.method = "average")
  ry <- rank(contribution, ties.method = "average")

  w <- weight
  w_sum <- sum(w, na.rm = TRUE)
  if (is.na(w_sum) || w_sum == 0) return(0)

  mx <- sum(w * rx, na.rm = TRUE) / w_sum
  my <- sum(w * ry, na.rm = TRUE) / w_sum

  num <- sum(w * (rx - mx) * (ry - my), na.rm = TRUE)
  den <- sqrt(sum(w * (rx - mx)^2, na.rm = TRUE) * sum(w * (ry - my)^2, na.rm = TRUE))

  if (is.na(den) || den == 0) return(0)
  rho <- num / den
  if (is.na(rho)) return(0)
  return(rho)
}


#' Evaluate Explainability Overfitting (Train vs Test SHAP)
#'
#' Compares the feature importance rankings between training and test sets
#' using Spearman and Kendall correlations. High correlation indicates that
#' the model's explanations generalize well (no overfitting of explanations).
#'
#' @param shap_raw_train Data.frame of SHAP values from training set.
#' @param shap_raw_test Data.frame of SHAP values from test set.
#'
#' @return A data.frame with per-model correlation metrics.
#' @export
evaluateExplainabilityOverfitting <- function(shap_raw_train, shap_raw_test) {
  models <- unique(shap_raw_train$model)

  purrr::map_dfr(models, function(mod) {
    train_sub <- shap_raw_train[shap_raw_train$model == mod, ]
    test_sub  <- shap_raw_test[shap_raw_test$model == mod, ]

    train_imp <- train_sub %>%
      dplyr::group_by(variable) %>%
      dplyr::summarize(imp_train = mean(abs(contribution), na.rm = TRUE), .groups = "drop")

    test_imp <- test_sub %>%
      dplyr::group_by(variable) %>%
      dplyr::summarize(imp_test = mean(abs(contribution), na.rm = TRUE), .groups = "drop")

    merged <- dplyr::inner_join(train_imp, test_imp, by = "variable")

    if (nrow(merged) < 3) {
      return(data.frame(model = mod, spearman_train_test = NA_real_, kendall_train_test = NA_real_))
    }

    spearman_rho <- suppressWarnings(stats::cor(merged$imp_train, merged$imp_test, method = "spearman"))
    kendall_tau  <- suppressWarnings(stats::cor(merged$imp_train, merged$imp_test, method = "kendall"))

    data.frame(
      model = mod,
      spearman_train_test = spearman_rho,
      kendall_train_test = kendall_tau
    )
  })
}

#' Evaluate Inter-Model Consistency (Jaccard and Spearman)
#'
#' Calculates the Jaccard similarity (on the top features) and the Spearman
#' correlation (on the global ranking) between different models. High agreement
#' indicates that the biological signature is robust and model-agnostic.
#'
#' @param shap_raw A data.frame of aggregated SHAP values (ideally from the test set).
#' @param top_n_jaccard Integer. Number of top features to use for Jaccard calculation.
#'
#' @return A data.frame with pairwise model comparisons, or NULL if < 2 models.
#' @export
evaluateInterModelConsistency <- function(shap_raw, top_n_jaccard = 20) {
  models <- unique(shap_raw$model)

  if (length(models) < 2) {
    message("  Only one model passed the cut-offs. Inter-model consistency skipped.")
    return(NULL)
  }

  model_signatures <- purrr::map(stats::setNames(models, models), function(mod) {
    df_mod <- shap_raw[shap_raw$model == mod, ]

    ranked <- df_mod %>%
      dplyr::group_by(variable) %>%
      dplyr::summarize(importance = mean(abs(contribution), na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(importance))

    list(
      model = mod,
      top_n = utils::head(ranked$variable, top_n_jaccard),
      scores = stats::setNames(ranked$importance, ranked$variable)
    )
  })

  pairs <- utils::combn(models, 2, simplify = FALSE)

  results <- purrr::map_dfr(pairs, function(p) {
    mod1 <- p[1]
    mod2 <- p[2]

    sig1 <- model_signatures[[mod1]]
    sig2 <- model_signatures[[mod2]]

    # Jaccard of Top N
    intersect_len <- length(intersect(sig1$top_n, sig2$top_n))
    union_len <- length(union(sig1$top_n, sig2$top_n))
    jaccard <- intersect_len / union_len

    # Spearman of full ranking
    common_vars <- intersect(names(sig1$scores), names(sig2$scores))
    rho <- if (length(common_vars) > 3) {
      suppressWarnings(stats::cor(sig1$scores[common_vars], sig2$scores[common_vars], method = "spearman"))
    } else {
      NA_real_
    }

    data.frame(
      Model_A = mod1,
      Model_B = mod2,
      Jaccard = jaccard,
      Spearman = rho
    )
  })

  return(results)
}
