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
#'   explainers <- createExplainers(final_models, df, "Group")
#'
#'   # Check class
#'   class(explainers$RF)
#'
createExplainers <- function(final_models, train_data, target_var) {
  message("Creating DALEX explainers...")
  positive_class <- levels(train_data[[target_var]])[2L]

  # Convert target to numeric 0/1 for DALEX
  y_numeric <- ifelse(
    train_data[[target_var]] == positive_class, 1L, 0L
  )

  # Remove target from features
  X <- train_data[, !colnames(train_data) %in% target_var, drop = FALSE]

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
#' @param target_var Target variable name.
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
                                     train_data,
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
  X_train <- train_data[, !colnames(train_data) %in% target_var, drop = FALSE]

  msg("Computing SHAP values...")

  shap <- run_map(explainers, function(exp) {
    shap_val <- DALEX::predict_parts(
      exp, new_observation = X_train,
      type = "shap", B = repetitions
    )
    shap_val$model <- exp$label
    shap_val
  })

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

  ## Top SHAP
  top_shap <- shap %>%
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

  ## Global consensus
  global_importance <- shap %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(
      mean_shap = mean(abs(contribution)),
      direction = mean(sign(contribution)),
      n_models  = dplyr::n_distinct(model),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(mean_shap))

  msg("Feature importance analysis complete.")

  list(
    permutation_raw = perm,
    permutation_top = top_perm,
    shap_raw = shap,
    shap_top = top_shap,
    global_importance = global_importance
  )
}
