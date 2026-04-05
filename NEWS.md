# RUMBLE 0.99.0

## Initial Bioconductor submission

### New features

* Complete consensus machine learning pipeline for microbiome biomarker
  discovery (`runPipeline()`).
* Accepts both `phyloseq` objects and count matrices as input
  (`prepareData()`).
* Tree-based consensus classification with Random Forest and XGBoost via
  `tidymodels`.
* Classical model-agnostic interpretability through SHAP values and
  permutation importance via `DALEX` (`computeFeatureImportance()`).
* Additional TreeSHAP interaction ecology stage with synergy,
  redundancy, and independence decomposition.
* Expanded publication-ready visualizations: metrics comparison, ROC
  curves, confusion matrices, SHAP consensus, SHAP distribution,
  permutation importance heatmap, S-R-I stacked decomposition, and
  ecological synergy/redundancy network assets.
* Optional removal of unclassified taxa (`remove_unclassified`
  parameter, default `FALSE`).
* Prevalence-based filtering and CLR transformation for compositional
  data.
