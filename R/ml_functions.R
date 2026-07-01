#' Internal Machine Learning Functions for RUMBLE
#'
#' These functions handle data splitting, recipe building, workflow
#' construction, hyperparameter tuning, and model evaluation. They are
#' called internally by \code{\link{RUMBLE}} and are not exported.
#'
#' @name ml-internal
#' @keywords internal
#' @return No return value. This object is internal and serves only for organization.
#'
#' @importFrom recipes recipe step_zv step_normalize step_corr
#'      all_predictors
# REMOVIDO: themis::step_downsample (usando Case Weights para balanceamento)
#' @importFrom parsnip rand_forest boost_tree logistic_reg
#'      nearest_neighbor set_engine set_mode fit
#' @importFrom workflows workflow add_model add_recipe
#' @importFrom tune tune tune_grid select_best finalize_workflow
#'      control_grid
#' @importFrom yardstick metric_set f_meas roc_auc accuracy
#'      bal_accuracy mcc
#' @importFrom rsample initial_split training testing vfold_cv
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom purrr map imap
#' @importFrom rlang sym
#' @importFrom dplyr bind_cols select bind_rows
#' @importFrom stats predict as.formula
#' @importFrom withr with_seed
#' @importFrom magrittr "%>%"
NULL


## ------------------------------------------------------------------
## Create Outer Cross-Validation Folds
## ------------------------------------------------------------------
#' @noRd
.createOuterFolds <- function(data, target_var, folds = 5L, seed = 42L) {
  if (!is.factor(data[[target_var]])) {
    data[[target_var]] <- as.factor(data[[target_var]])
  }

  withr::with_seed(seed, {
    folds_obj <- rsample::vfold_cv(
      data, v = folds,
      strata = !!rlang::sym(target_var)
    )
  })

  folds_obj
}

## ------------------------------------------------------------------
## Build preprocessing recipe
## ------------------------------------------------------------------
#' Build a preprocessing recipe for the ML pipeline.
#'
#' The recipe applies:
#' 1. Zero-variance filter (step_zv) to remove constant features.
#' 2. Centering and scaling (step_normalize) for KNN and ENET compatibility.
#'    Note: CLR-transformed data is already on a log-ratio scale, but
#'    step_normalize ensures all features have unit variance, which is
#'    critical for distance-based (KNN) and regularized (ENET) models.
#'    For tree-based models (RF, XGB), normalization is mathematically
#'    irrelevant but does not harm performance.
#' 3. Optional correlation filter (step_corr) if threshold is provided.
#'
#' \strong{Note on Class Imbalance:} Class imbalance is now handled via
#' cost-sensitive learning (case weights) rather than downsampling. This
#' preserves 100% of the original data, which is critical for small cohorts
#' (N < 200) where discarding samples destroys variance structure and reduces
#' statistical power.
#'
#' @noRd
.buildRecipe <- function(train_data, target_var,
                         balance_classes = TRUE, seed = 42L,
                         correlation_threshold = NULL) {

  rec <- recipes::recipe(
    stats::as.formula(paste(target_var, "~ .")),
    data = train_data
  )

  # Não precisamos mais do update_role para case_weight

  rec <- rec %>%
    recipes::step_zv(recipes::all_predictors()) %>%
    recipes::step_normalize(recipes::all_predictors())

  # Aplicando o Downsampling logo após a normalização
  if (balance_classes) {
    rec <- rec %>%
      themis::step_downsample(!!rlang::sym(target_var), seed = seed)
  }

  if (!is.null(correlation_threshold)) {
    rec <- rec %>%
      recipes::step_corr(recipes::all_predictors(), threshold = correlation_threshold)
  }

  return(rec)
}

## ------------------------------------------------------------------
## Define model workflows
## ------------------------------------------------------------------
#' @noRd
.buildWorkflows <- function(recipe, xgb_trees = 1000L, rf_trees = 500L, seed = 42L) {
  models <- list(
    RF = parsnip::rand_forest(
      mtry = tune::tune(),
      trees = rf_trees,
      min_n = tune::tune()
    ) %>%
      parsnip::set_engine("ranger", num.threads = 1, seed = seed) %>%
      parsnip::set_mode("classification"),

    XGB = parsnip::boost_tree(
      trees = xgb_trees,
      tree_depth = tune::tune(),
      learn_rate = tune::tune(),
      loss_reduction = tune::tune(),
      min_n = tune::tune()
    ) %>%
      parsnip::set_engine("xgboost", nthread = 1) %>%
      parsnip::set_mode("classification"),

    ENET = parsnip::logistic_reg(
      penalty = tune::tune(),
      mixture = tune::tune()
    ) %>%
      parsnip::set_engine("glmnet") %>%
      parsnip::set_mode("classification"),

    KNN = parsnip::nearest_neighbor(
      neighbors = tune::tune(),
      weight_func = tune::tune()
    ) %>%
      parsnip::set_engine("kknn") %>%
      parsnip::set_mode("classification")
  )

  purrr::map(models, function(mod) {
    workflows::workflow() %>%
      workflows::add_model(mod) %>%
      workflows::add_recipe(recipe)
  })
}

## ------------------------------------------------------------------
## Tune hyperparameters
## ------------------------------------------------------------------
#' @noRd
.tuneModels <- function(workflows, train_data, target_var,
                        grid_size = 10L, cv_folds = 5L,
                        n_cores = 1L, seed = 42L) {

  withr::with_seed(seed, {
    folds <- rsample::vfold_cv(
      train_data, v = cv_folds,
      strata = !!rlang::sym(target_var)
    )
  })

  metrics <- yardstick::metric_set(
    yardstick::f_meas,
    yardstick::roc_auc,
    yardstick::accuracy,
    yardstick::mcc
  )

  # Setup parallel backend for tuning if requested
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    message("Parallel tuning active: ", n_cores, " cores")
  }

  purrr::imap(workflows, function(wf, name) {
    message("  Tuning model: ", name)

    tryCatch({
      res <- tune::tune_grid(
        wf,
        resamples = folds,
        grid = grid_size,
        metrics = metrics,
        control = tune::control_grid(
          save_pred = TRUE,
          verbose = FALSE,
          allow_par = (n_cores > 1L)
        )
      )

      best_param <- tune::select_best(res, metric = "mcc")

      if (is.null(best_param) || nrow(best_param) == 0) {
        warning("Could not select best parameters for ", name)
        return(NULL)
      }

      final_wf <- tune::finalize_workflow(wf, best_param)

      list(
        tuning_results = res,
        best_params = best_param,
        final_workflow = final_wf
      )
    }, error = function(e) {
      warning("Tuning failed for ", name, ": ", e$message)
      return(NULL)
    })
  })
}


## ------------------------------------------------------------------
## Evaluate final models on test set
## ------------------------------------------------------------------
#' @noRd
.evaluateModels <- function(tuned_results, train_data,
                            test_data, target_var) {

  message("Evaluating models on test set...")
  metrics_fn <- yardstick::metric_set(
    yardstick::f_meas,
    yardstick::accuracy,
    yardstick::bal_accuracy,
    yardstick::mcc
  )

  ## Ensure consistent factor levels
  train_data[[target_var]] <- factor(train_data[[target_var]])
  test_data[[target_var]]  <- factor(
    test_data[[target_var]],
    levels = levels(train_data[[target_var]])
  )

  # Filter out failed models (NULLs from tuning)
  tuned_results <- tuned_results[!vapply(tuned_results, is.null, logical(1L))]

  purrr::imap(tuned_results, function(obj, name) {
    message("  Fitting final model: ", name)

    fitted_wf <- parsnip::fit(obj$final_workflow,
                              data = train_data)

    preds <- dplyr::bind_cols(
      sample_id = rownames(test_data), # <-- Nova linha garantindo o ID da amostra
      stats::predict(fitted_wf, test_data, type = "prob"),
      stats::predict(fitted_wf, test_data, type = "class"),
      test_data[, target_var, drop = FALSE]
    )

    perf <- metrics_fn(
      preds,
      truth = !!rlang::sym(target_var),
      estimate = .pred_class
    )

    # Calculate ROC AUC if probability column exists
    pos_class <- levels(train_data[[target_var]])[2L]
    col_prob  <- paste0(".pred_", pos_class)

    if (col_prob %in% colnames(preds)) {
      auc <- yardstick::roc_auc(
        preds,
        truth = !!rlang::sym(target_var),
        !!rlang::sym(col_prob),
        event_level = "second"
      )
      perf <- dplyr::bind_rows(perf, auc)
    }

    list(
      model_fit = fitted_wf,
      predictions = preds,
      metrics = perf
    )
  })
}
