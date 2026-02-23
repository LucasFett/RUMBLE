# RUMBLE 0.99.0

## Initial Bioconductor submission

### New features

* Complete consensus machine learning pipeline for microbiome biomarker
  discovery (`runPipeline()`).
* Accepts both `phyloseq` objects and count matrices as input
  (`prepareData()`).
* Four classification algorithms: Random Forest, XGBoost, Elastic Net,
  and KNN, all via `tidymodels`.
* Model-agnostic interpretability through SHAP values and permutation
  importance via `DALEX` (`computeFeatureImportance()`).
* Seven publication-ready visualizations: metrics comparison, ROC
  curves, confusion matrices, SHAP consensus, SHAP distribution,
  permutation importance barplot, and permutation importance heatmap.
* Optional removal of unclassified taxa (`remove_unclassified`
  parameter, default `FALSE`).
* Prevalence-based filtering and CLR transformation for compositional
  data.
