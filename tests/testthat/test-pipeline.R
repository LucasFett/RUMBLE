test_that("RUMBLE runs end-to-end and returns the new tree ecology layer", {
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

  # 2. Run the pipeline in minimalist mode
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
    verbose = FALSE
  )

  # 3. Verifications
  expect_type(results, "list")
  expect_true("models" %in% names(results))
  expect_true("plots" %in% names(results))
  expect_true("tree_ecology" %in% names(results))
  expect_true(all(names(results$models) %in% c("RF", "XGB")))
  expect_false(any(names(results$models) %in% c("ENET", "KNN")))
  expect_s3_class(results$plots$roc, "ggplot")
  expect_s3_class(results$plots$sri_decomposition, "ggplot")
  expect_s3_class(results$plots$sri_ecological_profile, "ggplot")
  expect_true(nrow(results$metrics) > 0)
  expect_true(is.data.frame(results$tree_ecology$profile))
  expect_true(is.data.frame(results$tree_ecology$synergy_edges))
  expect_true(is.data.frame(results$tree_ecology$redundancy_edges))
})
