test_that("RUMBLE pipeline returns consensus and model-specific interpretability outputs", {
  set.seed(123)
  n_samples <- 60
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
    tax_level = NULL,
    project_name = "Test_Run",
    output_dir = tempdir(),
    n_cores = 1,
    grid_size = 1,
    cv_folds = 2,
    shap_reps = 1,
    top_n = 2,
    verbose = FALSE,
    shap_model = "all"
  )

  expect_type(results, "list")
  expect_true(all(c("models", "plots", "importance", "selected_models") %in% names(results)))
  expect_true(nrow(results$metrics) > 0)

  expect_s3_class(results$plots$roc, "ggplot")
  expect_s3_class(results$plots$shap_spearman, "ggplot")
  expect_s3_class(results$plots$shap_mean, "ggplot")
  expect_s3_class(results$plots$shap_beeswarm, "ggplot")
  expect_true(is.list(results$plots$model_specific))

  expect_true(length(results$selected_models) >= 1)
  expect_identical(sort(names(results$plots$model_specific$shap_spearman)), sort(results$selected_models))
  expect_identical(sort(names(results$plots$model_specific$shap_mean)), sort(results$selected_models))
  expect_identical(sort(names(results$plots$model_specific$shap_beeswarm)), sort(results$selected_models))

  first_model <- results$selected_models[[1]]
  expect_s3_class(results$plots$model_specific$shap_spearman[[first_model]], "ggplot")
  expect_s3_class(results$plots$model_specific$shap_mean[[first_model]], "ggplot")
  expect_s3_class(results$plots$model_specific$shap_beeswarm[[first_model]], "ggplot")

  expect_true(all(c("spearman", "mean_shap") %in% names(results$importance$shap_top)))
  expect_true("model" %in% colnames(results$importance$shap_top$spearman))
  expect_true(all(results$selected_models %in% unique(results$importance$shap_top$spearman$model)))
})

test_that("shap_model = 'consensus' suppresses model-specific plots", {
  set.seed(456)
  n_samples <- 40
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
    tax_level = NULL,
    n_cores = 1,
    grid_size = 1,
    cv_folds = 2,
    shap_reps = 1,
    top_n = 2,
    verbose = FALSE,
    shap_model = "consensus"
  )

  expect_length(results$plots$model_specific$shap_spearman, 0)
  expect_length(results$plots$model_specific$shap_mean, 0)
  expect_length(results$plots$model_specific$shap_beeswarm, 0)
})

test_that("invalid shap_model names raise a clear error", {
  set.seed(789)
  n_samples <- 60
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
      tax_level = NULL,
      n_cores = 1,
      grid_size = 1,
      cv_folds = 2,
      shap_reps = 1,
      top_n = 2,
      verbose = FALSE,
      shap_model = "NOT_A_MODEL"
    ),
    "Invalid 'shap_model' selection"
  )
})
