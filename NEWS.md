## RUMBLE 3.1.0

### Major Methodological Corrections & New Features
* **Nested Cross-Validation Metric Consistency:** Fixed an optimization reporting bias in the Inner CV loop. The pipeline now rigorously extracts all performance metrics (AUC, MCC, etc.) exclusively from the exact hyperparameter configuration that won the fold. Previously, metrics could be artificially inflated due to independent `max()` aggregation across different `.config` sets.
* **Generalization Gap Accuracy:** With the refactored `summarizeInnerCV()` and updated extraction in `run_analysis.R`, the reported Generalization Gap now reflects perfect mathematical alignment between the exact model tuned in the inner loop and its performance on the out-of-fold test set.

## RUMBLE 3.0.1

* **Minor updates and fixes** for compatibility and robustness.

## RUMBLE 2.0.0

* **Mandatory `shap_data` Parameter:** The `shap_data` parameter is now mandatory and must be explicitly specified by the user. Previously, it defaulted to `"all"`. This change enforces conscious methodological decision-making: users must now choose between `"train"` (recommended for biomarker discovery), `"test"` (for prediction-focused analysis), or `"all"` (for maximum statistical power). This ensures greater methodological rigor and transparency in analytical objectives.

## RUMBLE 1.0.0

### Major Changes & New Features
* **Integrated Differential Abundance (DA):** Added `run_da` and `da_method` parameters to the main `RUMBLE()` interface. The pipeline now validates complex machine learning predictive patterns with classical DA methods (Wilcoxon, ANCOM-BC, or ALDEx2).
* **Total Interpretability Mode:** Changed the default value of `shap_data` to `"all"` (previously `"test"`). This maximizes statistical power for SHAP calculations and ensures perfect biological alignment between SHAP importance, group prevalence, and DA Log2FC.
* **Optimized Defaults:** Adjusted default `shap_reps` to `25` (down from 100), providing an optimal out-of-the-box balance between feature ranking stability and computation time.

## RUMBLE 0.99.0

### Reviewer-oriented improvements
* Added `shap_model` to the main `RUMBLE()` interface so reviewers and users can inspect SHAP results either as global consensus or isolated by selected model.
* Added model-specific SHAP, beeswarm, prevalence, and integrated biomarker plots for every model retained after `metric_cutoffs`.
* Added `selected_models` to the return object to make the interpretability inclusion rule explicit and auditable.
* Extended file export so `output_dir` now includes both consensus-level and model-specific SHAP artifacts.
* Improved argument documentation for `shap_data`, `shap_model`, and the updated return structure.
* Updated plotting helper documentation to support per-model interpretability views.
* Expanded automated tests to cover model-specific SHAP outputs, consensus-only mode, and invalid `shap_model` values.
* Rewrote the `README.md` to make installation, quick validation, reviewer inspection, and expected outputs clearer.

## RUMBLE 0.1.0

### New features
* Complete consensus machine learning pipeline for microbiome biomarker discovery (`RUMBLE()`).
* Accepts both `phyloseq` objects and count matrices as input (`prepareData()`).
* Four classification algorithms: Random Forest, XGBoost, Elastic Net, and KNN, all via `tidymodels`.
* Model-agnostic interpretability through SHAP values and permutation importance via `DALEX` (`computeFeatureImportance()`).
* Publication-ready visualizations for metrics, ROC curves, confusion matrices, SHAP summaries, integrated biomarker interpretation, and permutation importance.
* Optional removal of unclassified taxa (`remove_unclassified` parameter, default `FALSE`).
* Prevalence-based filtering and CLR transformation for compositional data.
