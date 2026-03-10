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
#'
#' @param explainers List of DALEX explainers (output of \code{createExplainers}).
#' @param train_data Training data.frame.
#' @param target_var Name of the target variable.
#' @param class_of_interest Character. The class of interest (must be the second
#'   factor level). SHAP values will be calculated with respect to this class.
#' @param top_n Number of top features to retain (default 20).
#' @param repetitions Number of repetitions (B) for importance calculation
#'   (default 10).
#' @param n_cores Number of cores for parallel processing (default 1).
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
#'      n_distinct bind_rows
#' @importFrom furrr future_map_dfr furrr_options
#' @importFrom future plan multisession sequential
#' @importFrom purrr map_dfr
#' @export
#'
#' @examples
#'   # ... (Assume 'explainers' and 'df' created as in createExplainers example)
#'   # For demonstration, we run with low repetitions
#'
#'   if (exists("explainers") && exists("df")) {
#'     importance <- computeFeatureImportance(
#'       explainers = explainers,
#'       train_data = df,
#'       target_var = "Group",
#'       top_n = 5,
#'       repetitions = 2,  # Low for speed
#'       n_cores = 1
#'     )
#'
#'     head(importance$global_importance)
#'   }
#'
computeFeatureImportance <- function(explainers,
                                     data,
                                     target_var,
                                     top_n = 20L,
                                     repetitions = 10L,
                                     n_cores = 1L,
                                     verbose = TRUE) {

  msg <- function(...) {
    if (verbose) message(...)
  }

  msg("------------------------------------------------------")
  msg("Feature importance analysis (RUMBLE)")
  msg("Models: ", length(explainers),
      " | Top N: ", top_n,
      " | Repetitions: ", repetitions)
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

  # --- 1. PERMUTATION IMPORTANCE ---
  msg("Computing permutation importance...")

  perm <- run_map(explainers, function(exp) {
    parts <- DALEX::model_parts(
      exp, type = "variable_importance", B = repetitions
    )
    # Filter out internal DALEX rows
    parts <- parts[!parts$variable %in% c("_full_model_", "_baseline_"), ]
    parts$model <- exp$label
    parts
  })

  # --- 2. SHAP VALUES ---
  # Remove target variable from training data for SHAP
  X_data <- data[, !colnames(data) %in% target_var, drop = FALSE]

  msg("Computing SHAP values for all observations (Global SHAP)...")

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
  # --- CLEANUP ---
  if (n_cores > 1) {
    future::plan(future::sequential)
  }

  # --- 3. AGGREGATE RESULTS ---
  msg("Aggregating results...")

  ## Top permutation
  top_perm <- perm %>%
    dplyr::group_by(model, variable) %>%
    dplyr::summarize(
      mean_dropout_loss = mean(dropout_loss, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(model, dplyr::desc(mean_dropout_loss)) %>%
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
