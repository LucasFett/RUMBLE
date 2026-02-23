test_that("RUMBLE works with synthetic matrix data", {
  # 1. Create synthetic data
  set.seed(123)
  n_samples <- 60
  n_taxa <- 50

  # False count matrix
  otu <- matrix(sample(0:100, n_samples * n_taxa, replace = TRUE),
                nrow = n_samples, ncol = n_taxa)
  rownames(otu) <- paste0("Sample", 1:n_samples)
  colnames(otu) <- paste0("Taxa", 1:n_taxa)

  # False metadata
  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  # 2. Run the pipeline in "minimalist" mode
  results <- RUMBLE(
    input = otu,
    metadata = meta,
    outcome_var = "Group",
    tax_level = NULL,
    project_name = "Test_Run",
    output_dir = tempdir(),
    n_cores = 1,
    grid_size = 1,
    cv_folds = 2,
    shap_reps = 1,
    top_n = 2,
    verbose = FALSE
  )

  # 3. Verifications (Expectations)
  expect_type(results, "list")
  expect_true("models" %in% names(results))
  expect_true("plots" %in% names(results))
  expect_s3_class(results$plots$roc, "ggplot")
  expect_true(nrow(results$metrics) > 0)
})
