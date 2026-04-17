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

## Initial Bioconductor submission

### New features

* Complete consensus machine learning pipeline for microbiome biomarker discovery (`RUMBLE()`).
* Accepts both `phyloseq` objects and count matrices as input (`prepareData()`).
* Four classification algorithms: Random Forest, XGBoost, Elastic Net, and KNN, all via `tidymodels`.
* Model-agnostic interpretability through SHAP values and permutation importance via `DALEX` (`computeFeatureImportance()`).
* Publication-ready visualizations for metrics, ROC curves, confusion matrices, SHAP summaries, integrated biomarker interpretation, and permutation importance.
* Optional removal of unclassified taxa (`remove_unclassified` parameter, default `FALSE`).
* Prevalence-based filtering and CLR transformation for compositional data.
