test_that("RUMBLE pipeline returns consensus and model-specific interpretability outputs", {
  set.seed(123)
  n_samples <- 100
  n_taxa <- 50

  otu <- matrix(
    sample(0:100, n_samples * n_taxa, replace = TRUE),
    nrow = n_samples,
    ncol = n_taxa
  )
  rownames(otu) <- paste0("Sample", seq_len(n_samples))
  colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  results <- RUMBLE(
    input = otu,
    metadata = meta,
    outcome_var = "Group",
    class_of_interest = "Disease",
    run_da = FALSE,
    tax_level = NULL,
    n_cores = 1,
    grid_size = 1,
    cv_folds = 2,
    shap_reps = 1,
    top_n = 2,
    verbose = FALSE,
    shap_data = "all",
    shap_model = "all",
    xgb_trees = 5,
    rf_trees = 5
  )

  expect_type(results, "list")
  expect_true(all(c("models", "plots", "importance", "selected_models", "cv_metrics") %in% names(results)))
  expect_true(nrow(results$metrics) > 0)

  expect_s3_class(results$plots$roc, "ggplot")
  expect_s3_class(results$plots$shap_spearman, "ggplot")
  expect_s3_class(results$plots$shap_beeswarm, "ggplot")
  expect_true(is.list(results$plots$model_specific))

  expect_true(length(results$selected_models) >= 1)
  expect_identical(sort(names(results$plots$model_specific$shap_spearman)), sort(results$selected_models))
  expect_identical(sort(names(results$plots$model_specific$shap_beeswarm)), sort(results$selected_models))

  first_model <- results$selected_models[[1]]
  expect_s3_class(results$plots$model_specific$shap_spearman[[first_model]], "ggplot")
  expect_s3_class(results$plots$model_specific$shap_beeswarm[[first_model]], "ggplot")

  expect_true("spearman" %in% names(results$importance$shap_top))
  expect_true("model" %in% colnames(results$importance$shap_top$spearman))
  expect_true(all(results$selected_models %in% unique(results$importance$shap_top$spearman$model)))
})


test_that("shap_model = 'consensus' suppresses model-specific plots", {
  set.seed(456)
  n_samples <- 100
  n_taxa <- 20

  otu <- matrix(
    sample(0:100, n_samples * n_taxa, replace = TRUE),
    nrow = n_samples,
    ncol = n_taxa
  )
  rownames(otu) <- paste0("Sample", seq_len(n_samples))
  colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  results <- RUMBLE(
    input = otu,
    metadata = meta,
    outcome_var = "Group",
    class_of_interest = "Disease",
    run_da = FALSE,
    tax_level = NULL,
    n_cores = 1,
    grid_size = 1,
    cv_folds = 2,
    shap_reps = 1,
    top_n = 2,
    verbose = FALSE,
    shap_data = "train",
    shap_model = "consensus",
    xgb_trees = 5,
    rf_trees = 5
  )

  expect_length(results$plots$model_specific$shap_spearman, 0)
  expect_length(results$plots$model_specific$shap_beeswarm, 0)
})


test_that("invalid shap_model names raise a clear error", {
  set.seed(789)
  n_samples <- 100
  n_taxa <- 15

  otu <- matrix(
    sample(0:100, n_samples * n_taxa, replace = TRUE),
    nrow = n_samples,
    ncol = n_taxa
  )
  rownames(otu) <- paste0("Sample", seq_len(n_samples))
  colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  expect_error(
    RUMBLE(
      input = otu,
      metadata = meta,
      outcome_var = "Group",
      class_of_interest = "Disease",
      run_da = FALSE,
      tax_level = NULL,
      n_cores = 1,
      grid_size = 1,
      cv_folds = 2,
      shap_reps = 1,
      top_n = 2,
      verbose = FALSE,
      shap_data = "test",
      shap_model = "NOT_A_MODEL",
      xgb_trees = 5,
      rf_trees = 5
    ),
    "Invalid 'shap_model' selection"
  )
})


test_that("RUMBLE runs differential abundance correctly when requested", {
  # Proteção de CI/CD: Este teste só roda se o usuário/servidor tiver as dependências do Bioconductor
  testthat::skip_if_not_installed("ANCOMBC")
  testthat::skip_if_not_installed("TreeSummarizedExperiment")

  set.seed(999)
  n_samples <- 100
  n_taxa <- 20

  otu <- matrix(
    sample(0:100, n_samples * n_taxa, replace = TRUE),
    nrow = n_samples,
    ncol = n_taxa
  )
  rownames(otu) <- paste0("Sample", seq_len(n_samples))
  colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  results_da <- RUMBLE(
    input = otu,
    metadata = meta,
    outcome_var = "Group",
    class_of_interest = "Disease",
    run_da = TRUE,
    da_method = "ancombc",
    tax_level = NULL,
    n_cores = 1,
    grid_size = 1,
    cv_folds = 2,
    shap_reps = 1,
    top_n = 2,
    verbose = FALSE,
    shap_data = "all",
    xgb_trees = 5,
    rf_trees = 5
  )

  expect_false(is.null(results_da$da_results))
  expect_true(is.data.frame(results_da$da_results))
  expect_true(all(c("Taxon", "logFC", "p_val", "adj_p_val") %in% colnames(results_da$da_results)))

  expect_s3_class(results_da$plots$biomarker_integrated_spearman, "ggplot")
})

test_that("RUMBLE correctly extracts inner CV metrics for generalization gap evaluation", {
  set.seed(321)
  n_samples <- 100
  n_taxa <- 20

  otu <- matrix(
    sample(0:100, n_samples * n_taxa, replace = TRUE),
    nrow = n_samples,
    ncol = n_taxa
  )
  rownames(otu) <- paste0("Sample", seq_len(n_samples))
  colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  results_cv <- RUMBLE(
    input = otu,
    metadata = meta,
    outcome_var = "Group",
    class_of_interest = "Disease",
    run_da = FALSE,
    tax_level = NULL,
    n_cores = 1,
    grid_size = 1,
    cv_folds = 2,
    shap_reps = 1,
    top_n = 2,
    verbose = FALSE,
    shap_data = "test",
    shap_model = "consensus",
    xgb_trees = 5,
    rf_trees = 5
  )

  # cv_metrics deve existir
  expect_false(is.null(results_cv$cv_metrics))
  expect_true(is.data.frame(results_cv$cv_metrics))

  # Deve conter as colunas esperadas do novo pipeline rigoroso
  expect_true(all(c("Model", "Fold", "Metric", "CV_Mean", "CV_StdErr") %in% colnames(results_cv$cv_metrics)))

  # A tabela não deve estar vazia
  expect_true(nrow(results_cv$cv_metrics) > 0)
})
