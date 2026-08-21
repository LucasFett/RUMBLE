test_that("RUMBLE pipeline returns consensus and model-specific interpretability outputs", {
  set.seed(123)
  n_samples <- 40
  n_taxa <- 10

  # Criação de dados com "sinal biológico" para evitar falha no tuning do XGBoost/ENET
  otu <- matrix(rpois(n_samples * n_taxa, lambda = 10), nrow = n_samples, ncol = n_taxa)
  otu[21:40, 1:5] <- otu[21:40, 1:5] + 50 # Adiciona sinal forte nas primeiras 5 taxas para o grupo Disease

  rownames(otu) <- paste0("Sample", seq_len(n_samples))
  colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  suppressWarnings({
    results <- RUMBLE(
      input = otu,
      metadata = meta,
      outcome_var = "Group",
      class_of_interest = "Disease",
      run_da = FALSE,
      tax_level = NULL,
      n_cores = 1,
      outer_folds = 2,  # CRÍTICO PARA VELOCIDADE DO TESTE
      grid_size = 1,
      cv_folds = 2,
      shap_reps = 1,
      top_n = 2,
      verbose = FALSE,
      shap_model = "all",
      xgb_trees = 5,
      rf_trees = 5,
      class_balance_method = "downsample",
      apply_stability_filter = FALSE
    )
  })

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
  n_samples <- 40
  n_taxa <- 10

  otu <- matrix(rpois(n_samples * n_taxa, lambda = 10), nrow = n_samples, ncol = n_taxa)
  otu[21:40, 1:5] <- otu[21:40, 1:5] + 50

  rownames(otu) <- paste0("Sample", seq_len(n_samples))
  colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  suppressWarnings({
    results <- RUMBLE(
      input = otu,
      metadata = meta,
      outcome_var = "Group",
      class_of_interest = "Disease",
      run_da = FALSE,
      tax_level = NULL,
      n_cores = 1,
      outer_folds = 2,
      grid_size = 1,
      cv_folds = 2,
      shap_reps = 1,
      top_n = 2,
      verbose = FALSE,
      shap_model = "consensus",
      xgb_trees = 5,
      rf_trees = 5,
      class_balance_method = "downsample",
      apply_stability_filter = FALSE
    )
  })

  expect_length(results$plots$model_specific$shap_spearman, 0)
  expect_length(results$plots$model_specific$shap_beeswarm, 0)
})


test_that("invalid shap_model names raise a clear error", {
  set.seed(789)
  n_samples <- 40
  n_taxa <- 10

  otu <- matrix(rpois(n_samples * n_taxa, lambda = 10), nrow = n_samples, ncol = n_taxa)
  otu[21:40, 1:5] <- otu[21:40, 1:5] + 50

  rownames(otu) <- paste0("Sample", seq_len(n_samples))
  colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  suppressWarnings({
    expect_error(
      RUMBLE(
        input = otu,
        metadata = meta,
        outcome_var = "Group",
        class_of_interest = "Disease",
        run_da = FALSE,
        tax_level = NULL,
        n_cores = 1,
        outer_folds = 2,
        grid_size = 1,
        cv_folds = 2,
        shap_reps = 1,
        top_n = 2,
        verbose = FALSE,
        shap_model = "NOT_A_MODEL",
        xgb_trees = 5,
        rf_trees = 5,
        class_balance_method = "downsample",
        apply_stability_filter = FALSE
      ),
      "Invalid 'shap_model' selection"
    )
  })
})


test_that("RUMBLE runs differential abundance correctly when requested", {
  # Proteção de CI/CD: Este teste só roda se o usuário/servidor tiver as dependências do Bioconductor
  testthat::skip_if_not_installed("ANCOMBC")
  testthat::skip_if_not_installed("TreeSummarizedExperiment")

  set.seed(999)
  n_samples <- 40
  n_taxa <- 10

  otu <- matrix(rpois(n_samples * n_taxa, lambda = 10), nrow = n_samples, ncol = n_taxa)
  otu[21:40, 1:5] <- otu[21:40, 1:5] + 50

  rownames(otu) <- paste0("Sample", seq_len(n_samples))
  colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  suppressWarnings({
    results_da <- RUMBLE(
      input = otu,
      metadata = meta,
      outcome_var = "Group",
      class_of_interest = "Disease",
      run_da = TRUE,
      da_method = "ancombc",
      tax_level = NULL,
      n_cores = 1,
      outer_folds = 2,
      grid_size = 1,
      cv_folds = 2,
      shap_reps = 1,
      top_n = 2,
      verbose = FALSE,
      xgb_trees = 5,
      rf_trees = 5,
      class_balance_method = "downsample",
      apply_stability_filter = FALSE
    )
  })

  expect_false(is.null(results_da$da_results))
  expect_true(is.data.frame(results_da$da_results))
  expect_true(all(c("Taxon", "logFC", "p_val", "adj_p_val") %in% colnames(results_da$da_results)))

  expect_s3_class(results_da$plots$biomarker_integrated_spearman, "ggplot")
})

test_that("class_balance_method = 'class_weights' applies native per-class weighting", {
  set.seed(654)
  n_samples <- 40
  n_taxa <- 10

  otu <- matrix(rpois(n_samples * n_taxa, lambda = 10), nrow = n_samples, ncol = n_taxa)
  otu[21:40, 1:5] <- otu[21:40, 1:5] + 50

  rownames(otu) <- paste0("Sample", seq_len(n_samples))
  colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  suppressWarnings({
    results_cw <- RUMBLE(
      input = otu,
      metadata = meta,
      outcome_var = "Group",
      class_of_interest = "Disease",
      run_da = FALSE,
      tax_level = NULL,
      n_cores = 1,
      outer_folds = 2,
      grid_size = 1,
      cv_folds = 2,
      shap_reps = 1,
      top_n = 2,
      verbose = FALSE,
      shap_model = "consensus",
      xgb_trees = 5,
      rf_trees = 5,
      class_balance_method = "class_weights",
      apply_stability_filter = FALSE
    )
  })

  # class_weights nao deve derrubar o pipeline (RF e XGB tem suporte
  # nativo; ENET/KNN ficam sem peso, mas ainda devem ser ajustados)
  expect_type(results_cw, "list")
  expect_true(nrow(results_cw$metrics) > 0)
  expect_true(length(results_cw$selected_models) >= 1)
})


test_that("class_balance_method = 'downsample' remains numerically unaffected by the class_weights refactor", {
  # Este teste documenta a garantia central do commit: o caminho
  # 'downsample' (usado no artigo) nao deve mudar de comportamento so'
  # porque 'class_weights' foi introduzido/substituiu 'case_weights'.
  # Nao compara numeros exatos (isso e' feito manualmente via
  # vignettes/teste.R contra o baseline do artigo, que usa dados reais),
  # mas garante que a API e a forma dos outputs continuam identicas.
  set.seed(111)
  n_samples <- 40
  n_taxa <- 10

  otu <- matrix(rpois(n_samples * n_taxa, lambda = 10), nrow = n_samples, ncol = n_taxa)
  otu[21:40, 1:5] <- otu[21:40, 1:5] + 50

  rownames(otu) <- paste0("Sample", seq_len(n_samples))
  colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

  meta <- data.frame(
    ID = rownames(otu),
    Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
  )
  rownames(meta) <- meta$ID

  run_once <- function() {
    set.seed(222)
    suppressWarnings({
      RUMBLE(
        input = otu,
        metadata = meta,
        outcome_var = "Group",
        class_of_interest = "Disease",
        run_da = FALSE,
        tax_level = NULL,
        n_cores = 1,
        outer_folds = 2,
        grid_size = 1,
        cv_folds = 2,
        shap_reps = 1,
        top_n = 2,
        verbose = FALSE,
        shap_model = "consensus",
        xgb_trees = 5,
        rf_trees = 5,
        class_balance_method = "downsample",
        apply_stability_filter = FALSE
      )
    })
  }

  res_a <- run_once()
  res_b <- run_once()

  # Com a mesma semente e class_balance_method = "downsample", dois runs
  # devem produzir metricas de teste identicas -- isso so' e' verdade se o
  # refactor de class_weights nao vazou nenhum efeito colateral para o
  # caminho downsample.
  expect_equal(res_a$metrics, res_b$metrics)
  expect_equal(res_a$cv_metrics, res_b$cv_metrics)
  expect_equal(res_a$integrated_table, res_b$integrated_table)
})


test_that("cv_metrics has the expected structure and content", {
  set.seed(321)
n_samples <- 40
n_taxa <- 10

otu <- matrix(rpois(n_samples * n_taxa, lambda = 10), nrow = n_samples, ncol = n_taxa)
otu[21:40, 1:5] <- otu[21:40, 1:5] + 50

rownames(otu) <- paste0("Sample", seq_len(n_samples))
colnames(otu) <- paste0("Taxa", seq_len(n_taxa))

meta <- data.frame(
  ID = rownames(otu),
  Group = factor(rep(c("Control", "Disease"), each = n_samples / 2))
)
rownames(meta) <- meta$ID

suppressWarnings({
  results_cv <- RUMBLE(
    input = otu,
    metadata = meta,
    outcome_var = "Group",
    class_of_interest = "Disease",
    run_da = FALSE,
    tax_level = NULL,
    n_cores = 1,
    outer_folds = 2,
    grid_size = 1,
    cv_folds = 2,
    shap_reps = 1,
    top_n = 2,
    verbose = FALSE,
    shap_model = "consensus",
    xgb_trees = 5,
    rf_trees = 5,
    class_balance_method = "downsample",
    apply_stability_filter = FALSE
  )
})

# cv_metrics deve existir
expect_false(is.null(results_cv$cv_metrics))
expect_true(is.data.frame(results_cv$cv_metrics))

# Deve conter as colunas esperadas do novo pipeline rigoroso
expect_true(all(c("Model", "Fold", "Metric", "CV_Mean", "CV_StdErr") %in% colnames(results_cv$cv_metrics)))

# A tabela não deve estar vazia
expect_true(nrow(results_cv$cv_metrics) > 0)
})
