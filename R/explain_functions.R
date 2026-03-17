#' Create DALEX Explainers for Trained Models
#'
#' Wraps each fitted tidymodels workflow into a DALEX explainer object.
#'
#' @param final_models List of model results from the pipeline. Each element
#'   must contain a \code{model_fit} component (the fitted workflow).
#' @param train_data Training data.frame (features + target).
#' @param target_var Name of the target variable.
#'
#' @return A list of DALEX explainer objects, one for each model.
#'
#' @importFrom DALEX explain
#' @importFrom purrr imap
#' @importFrom stats predict
#' @export
#'
#' @examples
#'   library(parsnip)
#'   library(workflows)
#'   library(dplyr)
#'
#'   # 1. Create synthetic data
#'   set.seed(123)
#'   df <- data.frame(
#'     X1 = rnorm(50),
#'     X2 = rnorm(50),
#'     Group = factor(rep(c("A", "B"), each = 25))
#'   )
#'
#'   # 2. Train a simple model
#'   model <- rand_forest(mode = "classification") %>%
#'     set_engine("ranger") %>%
#'     fit(Group ~ ., data = df)
#'
#'   # 3. Mock the RUMBLE output structure
#'   # (createExplainers expects a list where each element has $model_fit)
#'   final_models <- list(
#'     RF = list(model_fit = model)
#'   )
#'
#'   # 4. Create Explainers
#'   explainers <- createExplainers(final_models, df, "Group",
#'                                  class_of_interest = "B")
#'
#'   # Check class
#'   class(explainers$RF)
#'
createExplainers <- function(final_models, data, target_var,
                             class_of_interest) {
  message("Creating DALEX explainers...")

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

  # Remove target from features
  X <- data[, !colnames(data) %in% target_var, drop = FALSE]

  purrr::imap(final_models, function(model_obj, model_name) {
    fitted_wf <- model_obj$model_fit

    DALEX::explain(
      model = fitted_wf,
      data = X,
      y = y_numeric,
      label = model_name,
      type = "classification",
      verbose = FALSE, # Keep DALEX silent during creation
      predict_function = function(m, newdata) {
        pred <- stats::predict(m, newdata, type = "prob")
        pred[[paste0(".pred_", positive_class)]]
      }
    )
  })
}


#' Compute Feature Importance via SHAP and Permutation (Parallelized)
#'
#' Calculates SHAP values and permutation importance using parallel processing.
#' Supports both exact SHAP (via DALEX) and fast approximation (via fastshap).
#'
#' @param explainers List of DALEX explainers (output of \code{createExplainers}).
#' @param data Data.frame for SHAP predictions (can be train or test).
#' @param target_var Name of the target variable.
#' @param class_of_interest Character. The class of interest (must be the second
#'   factor level). SHAP values will be calculated with respect to this class.
#' @param top_n Number of top features to retain (default 20).
#' @param repetitions Number of repetitions (B) for importance calculation
#'   (default 10).
#' @param n_cores Number of cores for parallel processing (default 1).
#' @param shap_method Character. Method for SHAP calculation. Options are
#'   \code{"exact"} (default, uses DALEX for rigorous SHAP values) or
#'   \code{"fast"} (uses fastshap package for faster approximation).
#'   The exact method is recommended for publication-quality results.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{permutation_raw}: Raw permutation importance data.
#'   \item \code{permutation_top}: Top features by permutation importance.
#'   \item \code{shap_raw}: Raw SHAP values.
#'   \item \code{shap_top}: Top features by SHAP values.
#'   \item \code{global_importance}: Consensus importance across models.
#' }
#'
#' @importFrom DALEX model_parts predict_parts
#' @importFrom dplyr group_by summarize arrange slice_head ungroup desc
#'      n_distinct bind_rows mutate
#' @importFrom furrr future_map_dfr furrr_options
#' @importFrom future plan multisession sequential
#' @importFrom purrr map_dfr imap
#' @importFrom tidyr pivot_longer
#' @importFrom stats cor
#' @export
#'
#' @examples
#'   # ... (Assume 'explainers' and 'df' created as in createExplainers example)
#'   # For demonstration, we run with low repetitions
#'
#'   if (exists("explainers") && exists("df")) {
#'     importance <- computeFeatureImportance(
#'       explainers = explainers,
#'       data = df,
#'       target_var = "Group",
#'       class_of_interest = "B",
#'       top_n = 5,
#'       repetitions = 2,  # Low for speed
#'       n_cores = 1,
#'       shap_method = "exact"
#'     )
#'
#'     head(importance$global_importance)
#'   }
#'
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
  if (n_cores > 1) {
    msg("Initializing parallel backend with ", n_cores, " cores...")
    future::plan(future::multisession, workers = n_cores)
  } else {
    future::plan(future::sequential)
  }

  # --- MAPPER DEFINITION (Bypass furrr if n_cores == 1) ---
  run_map <- function(x, fn, seed = TRUE) {
    if (n_cores > 1) {
      furrr::future_map_dfr(x, fn,
                            .options = furrr::furrr_options(seed = seed))
    } else {
      purrr::map_dfr(x, fn)
    }
  }

  # --- PREPARE DATA FOR PERMUTATION AND SHAP ---
  X_data <- data[, !colnames(data) %in% target_var, drop = FALSE]
  y_numeric <- ifelse(data[[target_var]] == class_of_interest, 1L, 0L)

  # --- 1. PERMUTATION IMPORTANCE ---
  msg("Computing permutation importance...")

  perm <- run_map(explainers, function(exp) {
    parts <- DALEX::model_parts(
      exp,
      data = X_data,
      y = y_numeric,
      type = "variable_importance",
      B = repetitions
    )
    # DO NOT filter out _full_model_ yet, we need it to calculate Delta loss
    parts$model <- exp$label
    parts
  })

  # --- 2. SHAP VALUES ---
  # X_data already prepared above for permutation
  # For SHAP, we use the same X_data

  if (shap_method == "exact") {
    msg("Computing SHAP values for all observations (Global SHAP - EXACT method)...")

    # WE INVERT THE LOOP: Sequential over models, Parallel over patients!
    shap <- purrr::map_dfr(explainers, function(exp) {
      msg(paste0("  -> Extracting SHAP for model: ", exp$label))

      # run_map (parallel if n_cores > 1) is now applied to the patients
      run_map(seq_len(nrow(X_data)), function(i) {

        # Extract a single observation
        single_obs <- X_data[i, , drop = FALSE]

        # Calculate local SHAP for this specific patient
        shap_val <- DALEX::predict_parts(
          exp,
          new_observation = single_obs,
          type = "shap",
          B = repetitions
        )

        # Add necessary metadata for RUMBLE
        shap_val$model <- exp$label
        shap_val$feature_value <- as.numeric(sub(".*= ", "", as.character(shap_val$variable)))
        shap_val$observation_id <- i # Patient traceability

        return(shap_val)
      })
    })

    # Clean the variable name by removing the value part (e.g., "Taxon_A = 2.5" becomes "Taxon_A")
    shap$variable <- sub(" =.*", "", as.character(shap$variable))

  } else if (shap_method == "fast") {
    # Fast SHAP using fastshap package (approximation)
    msg("Computing SHAP values using FAST approximation method...")

    # Check if fastshap is installed
    if (!requireNamespace("fastshap", quietly = TRUE)) {
      stop("The 'fastshap' package is required for shap_method='fast'. ",
           "Install it with: install.packages('fastshap')")
    }

    shap <- purrr::map_dfr(explainers, function(exp) {
      msg(paste0("  -> Computing fast SHAP for model: ", exp$label))

      # Extract model and data
      fitted_model <- exp$model
      X_pred <- X_data

      # Define prediction function for fastshap
      pred_fn <- function(object, newdata) {
        stats::predict(object, newdata, type = "prob")[[paste0(".pred_", class_of_interest)]]
      }

      # Compute SHAP values using fastshap
      shap_vals <- fastshap::explain(
        object = fitted_model,
        X = X_pred,
        pred_wrapper = pred_fn,
        nsim = repetitions
      )

      # Convert to RUMBLE format
      shap_long <- tidyr::pivot_longer(
        data.frame(observation_id = seq_len(nrow(shap_vals)), shap_vals),
        cols = -observation_id,
        names_to = "variable",
        values_to = "contribution"
      )

      # PERFORMANCE OPTIMIZATION: Use vectorized join instead of row-by-row loop
      # Pivot original data to long format for efficient joining
      X_long <- X_pred %>%
        dplyr::mutate(observation_id = dplyr::row_number()) %>%
        tidyr::pivot_longer(
          cols = -observation_id,
          names_to = "variable",
          values_to = "feature_value"
        )

      # Join with SHAP values (vectorized operation in C++ backend)
      shap_long <- shap_long %>%
        dplyr::left_join(X_long, by = c("observation_id", "variable")) %>%
        dplyr::mutate(model = exp$label)

      return(shap_long)
    })
  }

  # --- CLEANUP ---
  if (n_cores > 1) {
    future::plan(future::sequential)
  }

  # --- 3. AGGREGATE RESULTS ---
  msg("Aggregating results...")

  ## Top permutation
  top_perm <- perm %>%
    dplyr::group_by(model) %>%
    dplyr::mutate(
      # Calculate full model loss for each model
      full_model_loss = mean(
        dropout_loss[variable == "_full_model_"],
        na.rm = TRUE
      )
    ) %>%
    dplyr::ungroup() %>%
    # Remove internal DALEX rows AFTER calculating full_model_loss
    dplyr::filter(!variable %in% c("_full_model_", "_baseline_")) %>%
    dplyr::mutate(
      # Calculate Delta loss (relative importance)
      delta_loss = dropout_loss - full_model_loss
    ) %>%
    dplyr::group_by(model, variable) %>%
    dplyr::summarize(
      mean_delta_loss = mean(delta_loss, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(model, dplyr::desc(mean_delta_loss)) %>%
    dplyr::group_by(model) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()

  ## --- 3.1 Spearman Metric ---
  ## Top SHAP (Spearman)
  top_shap_spearman <- shap %>%
    dplyr::group_by(model, variable) %>%
    dplyr::summarize(
      mean_abs_contribution = mean(abs(contribution), na.rm = TRUE),
      direction = sign(stats::cor(feature_value, contribution, method = "spearman")),
      .groups = "drop"
    ) %>%
    dplyr::arrange(model, dplyr::desc(mean_abs_contribution)) %>%
    dplyr::group_by(model) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()

  ## Global consensus (Spearman)
  global_importance_spearman <- shap %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      mean_shap = mean(abs(contribution), na.rm = TRUE),
      direction = sign(stats::cor(feature_value, contribution, method = "spearman")),
      n_models  = dplyr::n_distinct(model),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_shap))

  ## --- 3.2 Mean SHAP Metric ---
  ## Top SHAP (Mean SHAP)
  top_shap_mean <- shap %>%
    dplyr::group_by(model, variable) %>%
    dplyr::summarize(
      mean_abs_contribution = mean(abs(contribution), na.rm = TRUE),
      direction = mean(sign(contribution), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(model, dplyr::desc(mean_abs_contribution)) %>%
    dplyr::group_by(model) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()

  ## Global consensus (Mean SHAP)
  global_importance_mean <- shap %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      mean_shap = mean(abs(contribution), na.rm = TRUE),
      direction = mean(sign(contribution), na.rm = TRUE),
      n_models  = dplyr::n_distinct(model),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_shap))

  msg("Feature importance analysis complete.")

  list(
    permutation_raw = perm,
    permutation_top = top_perm,
    shap_raw = shap,
    shap_top = list(
      spearman = top_shap_spearman,
      mean_shap = top_shap_mean
    ),
    global_importance = list(
      spearman = global_importance_spearman,
      mean_shap = global_importance_mean
    )
  )
}
