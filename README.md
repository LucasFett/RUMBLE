# RUMBLE

## Reproducible Understanding of Microbiome Biomarkers with Leakage-free Explainability

**RUMBLE** is an R package for microbiome biomarker discovery designed around **predictive robustness**, **model-agnostic interpretability**, and **analytical reproducibility**. The package integrates multiple supervised classification algorithms within the **tidymodels** ecosystem, applies compositional-data-aware preprocessing, and returns a structured set of outputs suitable for both exploratory analysis and formal technical review.

The central design principle of RUMBLE is the explicit separation between two analytical layers. The first layer is **predictive performance**, assessed on a held-out test set using multiple evaluation metrics. The second layer is **interpretability**, constructed through **SHAP** values and **permutation importance** only for models that are considered eligible for interpretation. This design allows users to impose objective quality filters before generating biomarker interpretations, thereby reducing the risk of explaining weak or unstable models.

In its current form, the package also includes a reviewer-oriented feature that is particularly important for methodological inspection: the ability to inspect **model-isolated results**. By default, `shap_model = "all"` generates, in addition to the consensus plots, a set of outputs specific to each model that passed `metric_cutoffs`. This allows reviewers to verify whether SHAP profiles are convergent across algorithms or whether the interpretative signature is strongly model-dependent.

---

## Contents

| Section | Description |
|---|---|
| Overview | Package goals and methodological design |
| Installation | Development installation and suggested dependencies |
| Quick start | Minimal reproducible example |
| Input and output structure | Expected data formats and return object |
| Interpretability | Global consensus, per-model profiles, and metric-based filtering |
| Nested CV & Generalization Gap | Inner CV metrics and validation of overfitting |
| Exported artifacts | Files written to `output_dir` |
| Reviewer checklist | How to validate the package quickly |
| Testing and development | Automated testing and package inspection |

---

## Overview

The main pipeline is executed through `RUMBLE()`. At a high level, the function prepares the data, splits the dataset into training and test partitions, tunes and fits multiple models, evaluates out-of-sample performance, and then computes interpretability outputs using SHAP and permutation-based importance.

The four algorithms currently integrated are **Random Forest (`RF`)**, **XGBoost (`XGB`)**, **Elastic Net (`ENET`)**, and **K-Nearest Neighbors (`KNN`)**. The modeling layer relies on the **tidymodels** ecosystem [1], whereas interpretability is primarily built on **DALEX** [2]. An approximate SHAP route using **fastshap** is also supported when faster exploratory runs are desired [3].

The broader purpose of the package is to provide a unified interface for a workflow that is often implemented by researchers through fragmented, manually maintained scripts. RUMBLE organizes this process into a standardized return object containing metrics, importance tables, plots, and the intermediate data objects needed for analytical auditing.

---

## Installation

### Development version

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("LucasFett/RUMBLE")
```

### Suggested dependencies for full coverage

Some modeling engines are provided through suggested packages. To ensure that all algorithms and examples run with full coverage, the following installation is recommended:

```r
install.packages(c("tidymodels", "glmnet", "kknn", "ranger", "xgboost"))
```

### Bioconductor installation

Once the package becomes available through Bioconductor, installation can be performed as follows:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("RUMBLE")
```

---

## Quick start

The example below illustrates a minimal workflow using the `phyloseq` object distributed in `inst/extdata`.

```r
library(RUMBLE)
library(phyloseq)

ps_path <- system.file(
  "extdata",
  "PRJEB38465_phyloseq_com_metadados_completos.rds",
  package = "RUMBLE"
)

ps <- readRDS(ps_path)

set.seed(42)
results <- RUMBLE(
  input = ps,
  outcome_var = "ses",
  class_of_interest = "Low",
  tax_level = "Genus",
  output_dir = "rumble_results",
  verbose = TRUE
)

results$metrics
results$plots$roc
results$plots$shap_spearman
```

This default workflow produces performance metrics, ROC curves, confusion matrices, SHAP plots, permutation importance summaries, and integrated biomarker interpretations.

---

## Input structure

The package accepts two primary input formats.

| Input type | Description | Requirements |
|---|---|---|
| `phyloseq` | Full object containing abundances and sample metadata | The outcome variable must be present in `sample_data` |
| `matrix` / `data.frame` | Numeric matrix with samples in rows and taxa in columns | `metadata` must be supplied separately with compatible `rownames` |

The outcome must currently be **binary**. The `class_of_interest` argument is mandatory because it explicitly defines which class should be treated as the positive event for interpretation. This ensures that SHAP directionality is semantically anchored throughout the entire pipeline.

---

## Key arguments

The table below summarizes the most relevant arguments for practical use and reviewer inspection.

| Argument | Analytical role | Default |
|---|---|---|
| `tax_level` | Aggregates taxa at a chosen taxonomic rank | `NULL` |
| `class_balance_method` | Controls class imbalance handling (`"downsample"`, `"class_weights"`, or `"none"`) | `"downsample"` |
| `apply_stability_filter` | Filters unstable features from final interpretations based on `min_fold_frequency` | `FALSE` |
| `grid_size` | Controls tuning intensity | `30` |
| `shap_reps` | Number of repetitions for SHAP and permutation importance | `25` |
| `shap_method` | SHAP strategy (`"exact"` or `"fast"`) | `"exact"` |
| `metric_cutoffs` | Minimum quality filter for selecting interpretable models | `NULL` |
| `shap_model` | Controls which model-specific SHAP profiles are generated | `"all"` |
| `output_dir` | Directory used for automatic export of plots and tables | `NULL` |

**Note:** SHAP importance is always computed on both the training and test partitions internally; there is no `shap_data` argument to choose between them. The test-set signature is used as the primary consensus source when available (`results$importance$global_importance_test`), and the training-set signature is retained for the overfitting diagnostic (`results$shap_overfitting`).

---

## Main return object

`RUMBLE()` returns a structured list designed for inspection, programmatic reuse, and formal review.

| Element | Content |
|---|---|
| `models` | Final fitted models and test-set predictions |
| `metrics` | Summary table of model performance metrics across outer test folds |
| `cv_metrics` | Strict, configuration-matched inner loop CV metrics for assessing Generalization Gap |
| `importance` | Raw SHAP (train/test), model-specific SHAP summaries, weighted global consensus, and permutation importance |
| `plots` | Performance and interpretability plots (see below) |
| `selected_models` | Models retained for interpretability after applying `metric_cutoffs` |
| `hyperparameters` | Best hyperparameter configuration selected per model, per outer fold |
| `shap_overfitting` | Train-vs-test SHAP ranking correlation per model (explanation stability diagnostic) |
| `model_consistency` | Pairwise inter-model agreement (Jaccard similarity, Spearman correlation) on SHAP rankings |
| `da_results` | Differential abundance results (only if `run_da = TRUE`) |
| `integrated_table` | Final consensus biomarker table (SHAP magnitude, direction, prevalence, DA effect size) |
| `shap_consistency` | Inter-fold SHAP ranking stability per model |
| `feature_frequency` | Per-model and consensus top-N selection frequency across outer folds (always unfiltered, regardless of `apply_stability_filter`) |
| `hyperparameter_stability` | Variance of selected hyperparameters across outer folds |
| `metrics_per_fold` | Raw per-fold performance metrics (before summarization into `metrics`) |

In practice, this means that the interpretability stage does not need to include every trained model. The `selected_models` object explicitly records which algorithms were considered eligible for the final interpretation layer.

---

## Interpretability: global consensus and model-isolated profiles

Biomarker interpretation is organized at two complementary scales.

The first is the **global consensus** scale, which is useful when the goal is to communicate an aggregated microbiome signature. At this level, the package returns plots such as `results$plots$shap_spearman` and `results$plots$biomarker_integrated_spearman`.

The second is the **model-isolated** scale, which is essential for critical review. At this level, the package returns named lists inside `results$plots$model_specific`, making it possible to compare the behavior of `RF`, `XGB`, `ENET`, and `KNN` whenever those models pass the quality filter.

### Behavior of `shap_model`

| `shap_model` value | Effect |
|---|---|
| `"all"` | Generates isolated outputs for every model selected for SHAP interpretation |
| `"consensus"` | Generates only the consensus-level plots |
| `c("RF", "XGB")` | Generates isolated outputs only for the explicitly requested models |

### Example 1: generate consensus and model-specific profiles

```r
results <- RUMBLE(
  input = ps,
  outcome_var = "ses",
  class_of_interest = "Low",
  tax_level = "Genus",
  metric_cutoffs = c(roc_auc = 0.75),
  shap_model = "all"
)

names(results$plots$model_specific$shap_spearman)
results$plots$model_specific$shap_spearman$RF
results$plots$model_specific$shap_beeswarm$XGB
```

### Example 2: consensus-only inspection

```r
results <- RUMBLE(
  input = ps,
  outcome_var = "ses",
  class_of_interest = "Low",
  tax_level = "Genus",
  shap_model = "consensus"
)
```

### Example 3: restrict inspection to specific models

```r
results <- RUMBLE(
  input = ps,
  outcome_var = "ses",
  class_of_interest = "Low",
  tax_level = "Genus",
  metric_cutoffs = c(roc_auc = 0.70, f_meas = 0.60),
  shap_model = c("RF", "XGB")
)
```

---

## Quality filtering with `metric_cutoffs`

The `metric_cutoffs` argument determines which models are allowed to contribute to the interpretability stage. This is important because it separates **trained models** from **interpretable models**. In other words, all models are still evaluated in `results$metrics`, but only models that satisfy the user-defined thresholds are allowed to contribute to SHAP tables and plots.

This behavior is especially valuable for review, because it makes the inclusion rule behind the consensus explicit. If no model satisfies the specified thresholds, the pipeline stops with a clear message instead of silently producing incoherent interpretations.

```r
results <- RUMBLE(
  input = ps,
  outcome_var = "ses",
  class_of_interest = "Low",
  tax_level = "Genus",
  metric_cutoffs = c(roc_auc = 0.80, bal_accuracy = 0.70),
  shap_model = "all"
)

results$selected_models
```

---

## Exported artifacts in `output_dir`

When `output_dir` is defined, the package automatically writes plots and tables to disk. In the current version, export also includes model-specific artifacts, which makes the package easier to audit independently of the in-memory R object.

| Artifact type | Example files |
|---|---|
| Global metrics | `*_model_metrics.tsv` |
| Inner CV metrics | `*_cv_metrics.tsv` |
| Consensus SHAP | `*_shap_global_spearman.tsv` |
| Model-specific SHAP | `*_shap_model_spearman.tsv` |
| Permutation importance | `*_permutation_top.tsv` |
| Global figures | `*_roc.png`, `*_cm.png`, `*_oof_density.png`, `*_shap_spearman.png` |
| Model-specific figures | `*_RF_shap_spearman.png`, `*_XGB_shap_beeswarm.png`, and related outputs |

This organization makes the output directory substantially more reviewer-friendly because direct inspection no longer depends on manual extraction from the returned object.

---

## Nested Cross-Validation & Generalization Gap

RUMBLE 3.1.0 enforces strict mathematical alignment in the extraction of metrics from the **Nested Cross-Validation** inner loop. When the hyperparameter tuning step selects the winning configuration (e.g., the one that maximizes MCC), RUMBLE extracts all other metrics (AUC, Accuracy, F1, etc.) exclusively from that exact same configuration (`.config`).

This prevents "Optimization Reporting Bias" (where one might incorrectly report the maximum AUC from configuration A alongside the maximum MCC from configuration B). The resulting `results$cv_metrics` table provides an accurate estimation of inner-fold performance. Comparing these inner metrics with the outer-fold test metrics (`results$metrics`) allows researchers to rigorously assess the **Generalization Gap** and diagnose potential model overfitting before deriving biomarker interpretations.

---

## User checklist

A fast technical review can be carried out in only a few steps. First, run the minimal example with reduced `grid_size` and `shap_reps`. Then verify that `results$metrics` was generated and that `results$selected_models` reflects the thresholds imposed through `metric_cutoffs`. Next, compare the global plots against the model-specific plots stored in `results$plots$model_specific`. Finally, inspect the files written to `output_dir` and verify that both consensus-level and model-isolated artifacts are present.

| Review question | Where to inspect |
|---|---|
| Were the models trained and evaluated correctly? | `results$metrics`, `results$plots$roc`, `results$plots$cm` |
| Was the quality filter applied correctly? | `results$selected_models` |
| Does the SHAP consensus depend disproportionately on one algorithm? | Compare `results$plots$shap_*` with `results$plots$model_specific$*` |
| Are the exported outputs sufficient for external auditing? | Inspect `output_dir` |
| Does the package include user-facing documentation and tutorial material? | `README.md`, `vignettes/RUMBLE_tutorial.Rmd` |

---

## Practical guidance

For exploratory analyses, it is reasonable to reduce `grid_size` and `shap_reps` in order to speed up iteration. For manuscript-oriented analyses, external review, or stronger reproducibility requirements, the recommended approach is to return to the defaults or use even more conservative settings, while keeping `shap_method = "exact"` whenever computationally feasible.

It is also advisable not to interpret the consensus output as sufficient biological truth by itself. The model-specific layer was added precisely so that users can distinguish between **signals that are consistent across algorithms** and **signals that are strongly induced by a single modeling strategy**.

---

## Citation

If you use RUMBLE in your research, please cite:

> Fett, L. (2026). *RUMBLE: Reproducible Understanding of Microbiome Biomarkers with Leakage-free Explainability*. R package version 0.99.0.

---

## License

The package is distributed under the **GPL-3** license.

---

## References

[1]: https://www.tidymodels.org/ "tidymodels"
[2]: https://modeloriented.github.io/DALEX/ "DALEX"
[3]: https://bgreenwell.github.io/fastshap/ "fastshap"
