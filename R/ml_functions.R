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
#' @importFrom themis step_downsample
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
#' \strong{Note on Class Imbalance:} Class imbalance can be handled via
#' three methods controlled by \code{class_balance_method}:
#' \itemize{
#'   \item \code{"downsample"}: Reduces majority class to match minority class size.
#'   \item \code{"class_weights"}: Applies native per-class weighting at the
#'     engine level (RF via \code{ranger::class.weights}, XGB via
#'     \code{xgboost::scale_pos_weight}), preserving 100\% of the training
#'     data. ENET and KNN have no native per-class weighting mechanism in
#'     their respective engines and are fit unweighted under this mode; see
#'     \code{.buildWorkflows()} for details.
#'   \item \code{"none"}: No class balancing applied (useful for already-balanced data).
#' }
#'
#' @noRd
.buildRecipe <- function(train_data, target_var,
                         class_balance_method = c("downsample", "class_weights", "none"),
                         seed = 42L,
                         correlation_threshold = NULL) {

  class_balance_method <- match.arg(class_balance_method)

  # Criar a fórmula com todos os preditores. Ao contrário da abordagem
  # anterior baseada em case weights, "class_weights" não precisa de
  # nenhuma coluna extra nem de update_role() aqui -- os pesos por classe
  # são passados diretamente aos engines (ranger/xgboost) em
  # .buildWorkflows(), então a recipe permanece igual em todos os modos.
  rec <- recipes::recipe(
    stats::as.formula(paste(target_var, "~ .")),
    data = train_data
  )

  rec <- rec %>%
    recipes::step_zv(recipes::all_predictors()) %>%
    recipes::step_normalize(recipes::all_predictors())

  if (class_balance_method == "downsample") {
    rec <- rec %>%
      themis::step_downsample(!!rlang::sym(target_var), seed = seed)
  }
  # "class_weights": nada a fazer na recipe; os pesos por classe são
  # aplicados diretamente nos model spec's engine args, em .buildWorkflows().
  # "none": nada a fazer.

  if (!is.null(correlation_threshold)) {
    rec <- rec %>%
      recipes::step_corr(recipes::all_predictors(), threshold = correlation_threshold)
  }

  return(rec)
}

## ------------------------------------------------------------------
## Define model workflows
## ------------------------------------------------------------------
#' Build model workflows for the four base learners.
#'
#' \strong{Note on \code{class_balance_method = "class_weights"}:} weights are
#' applied natively at the engine level, not via the generic tidymodels
#' case-weights framework, because that framework (1) requires
#' \code{recipes}/\code{workflows} plumbing that proved fragile in practice
#' (e.g. \code{recipes::recipe()} does not accept \code{~ . - col} formula
#' syntax to exclude a weights column, and role-based workarounds add
#' complexity), and (2) does not have uniform support across the four engines
#' used here (notably \code{kknn}, which does not support case weights at
#' all). Native per-class weighting is simpler and avoids that whole class of
#' bug, at the cost of only being available for engines that expose a
#' class-level (not just observation-level) weighting argument:
#' \itemize{
#'   \item \strong{RF (ranger)}: \code{class.weights}, a numeric vector named
#'     by outcome factor level, applied via \code{set_engine()}.
#'   \item \strong{XGB (xgboost)}: \code{scale_pos_weight}, the
#'     negative/positive class count ratio, applied via \code{set_engine()}.
#'   \item \strong{ENET (glmnet)} and \strong{KNN (kknn)}: neither engine
#'     exposes a native per-class weighting argument (glmnet only supports
#'     per-observation weights, which would require the same fragile
#'     case-weights plumbing we are moving away from; kknn supports neither).
#'     These two models are fit \emph{without} imbalance correction under
#'     \code{"class_weights"} mode. A warning is emitted once per call to
#'     make this explicit rather than silent.
#' }
#'
#' @noRd
.buildWorkflows <- function(recipe, xgb_trees = 1000L, rf_trees = 500L, seed = 42L,
                            class_balance_method = "downsample",
                            class_weights_rf = NULL,
                            scale_pos_weight_xgb = NULL) {

  use_class_weights <- identical(class_balance_method, "class_weights")

  if (use_class_weights && (is.null(class_weights_rf) || is.null(scale_pos_weight_xgb))) {
    warning(
      "class_balance_method = 'class_weights' was requested but ",
      "class_weights_rf/scale_pos_weight_xgb were not supplied. ",
      "RF and XGB will be fit WITHOUT class weighting for this fold.",
      call. = FALSE
    )
    use_class_weights <- FALSE
  }

  if (use_class_weights) {
    message(
      "  class_balance_method = 'class_weights': applying native weighting ",
      "to RF (ranger::class.weights) and XGB (xgboost::scale_pos_weight). ",
      "ENET and KNN have no native per-class weighting mechanism and are ",
      "fit unweighted under this mode."
    )
  }

  rf_spec <- parsnip::rand_forest(
    mtry = tune::tune(),
    trees = rf_trees,
    min_n = tune::tune()
  ) %>%
    parsnip::set_mode("classification")

  rf_spec <- if (use_class_weights) {
    rf_spec %>% parsnip::set_engine(
      "ranger", num.threads = 1, seed = seed,
      class.weights = class_weights_rf
    )
  } else {
    rf_spec %>% parsnip::set_engine("ranger", num.threads = 1, seed = seed)
  }

  xgb_spec <- parsnip::boost_tree(
    trees = xgb_trees,
    tree_depth = tune::tune(),
    learn_rate = tune::tune(),
    loss_reduction = tune::tune(),
    min_n = tune::tune()
  ) %>%
    parsnip::set_mode("classification")

  xgb_spec <- if (use_class_weights) {
    xgb_spec %>% parsnip::set_engine(
      "xgboost", nthread = 1, scale_pos_weight = scale_pos_weight_xgb
    )
  } else {
    xgb_spec %>% parsnip::set_engine("xgboost", nthread = 1)
  }

  models <- list(
    RF = rf_spec,
    XGB = xgb_spec,

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

  purrr::imap(models, function(mod, model_name) {
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
