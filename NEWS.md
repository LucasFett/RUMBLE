## RUMBLE 3.2.0

### Major Methodological Corrections & New Features
* **Class Imbalance Handling Flexibility:** Introduced `class_balance_method` to `RUMBLE()`, allowing users to choose between `"downsample"` (default), `"class_weights"` (native per-class weighting at the engine level -- `ranger::class.weights` for RF and `xgboost::scale_pos_weight` for XGB -- preserving all training data; ENET and KNN have no native per-class weighting mechanism and are fit unweighted under this mode), or `"none"`. This ensures robust performance on imbalanced datasets without unnecessarily discarding data, especially for small cohorts.
* **Feature Stability Filtering:** Added `apply_stability_filter` to `RUMBLE()`. When enabled, only features appearing in the top-N for at least `min_fold_frequency` fraction of folds are included in the final integrated biomarker table and consensus plots. The raw `feature_frequency` table remains unfiltered for comprehensive auditing.
* **CLR / Prevalence-Filter Ordering Fix:** Corrected the intra-fold preprocessing order. CLR transformation is now applied strictly *after* the intra-fold prevalence filter, so the per-sample CLR geometric mean is computed only over taxa that survive filtering. Previously, CLR was computed before filtering, so taxa later discarded still contributed to the log-ratio denominator -- a substantive computational inconsistency with the documented anti-leakage design, not merely a documentation gap. All downstream CLR-transformed values (train and test, per outer fold) are affected; the `.apply_clr()` pseudocount is unchanged (`1e-6`, still hardcoded -- see Known Limitations below).
* **Weighted Direction Metric in Global Consensus:** The `SHAP_Rho` direction metric in the weighted global consensus (`global_importance_test`/`global_importance_train`) is now aggregated using the same MCC-derived per-model weight already applied to the magnitude metric (`mean_abs_contribution`). Previously, magnitude was weighted by model performance but direction was not, an inconsistency within the same consensus step.
* **Removed Unused Consensus Path:** Removed the `model_weights` parameter and the associated `global_importance` field from `computeFeatureImportance()`, along with the per-fold performance-weight calculation in `run_analysis.R` that fed it. Neither was consumed downstream -- the actual reported consensus has always been computed independently from raw per-fold SHAP values. The only remaining performance-weighting mechanism is the MCC-based `global_weights` used in the global consensus step, matching the manuscript's Eq. 5.
* **Documentation Debt -- `shap_data` Removal:** Backfilling documentation for an earlier, previously undocumented change: the `shap_data` parameter (`"train"`/`"test"`/`"all"`) no longer exists. `RUMBLE()` now always computes SHAP on both training and test partitions; the test-set signature is used as the primary consensus source when available (gold standard), and the training-set signature is retained for the train/test overfitting diagnostic (`results$shap_overfitting`). Examples referencing `shap_data` in older documentation are stale and have been corrected.

### Known Limitations (tracked for a future release)
* `pseudocount` in `.apply_clr()` remains hardcoded (`1e-6`) in three call sites rather than exposed as a `RUMBLE()` parameter.
* SHAP/permutation importance computation (`computeFeatureImportance()`) falls back to single-core execution for the *entire* model set whenever an XGBoost explainer is present, even though only the XGBoost path has a serialization constraint; RF/ENET/KNN could in principle still run in parallel.

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
