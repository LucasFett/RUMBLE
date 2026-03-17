# RUMBLE
## Robust Unified Microbiome Biomarker Learning Engine

RUMBLE is an R/Bioconductor package that implements a consensus machine learning framework for identifying microbial biomarkers from high-throughput sequencing data.

---

# Overview

RUMBLE integrates multiple classification algorithms (**Random Forest**, **XGBoost**, **Elastic Net**, and **KNN**) via the *tidymodels* ecosystem and provides model-agnostic interpretability through SHAP values and permutation importance (via *DALEX*).

The consensus approach ranks biomarkers across all models, using Spearman correlation to correctly map SHAP directionality, providing robust and reproducible feature selection for microbiome research.

---

# Installation

## Install from Bioconductor (when available)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RUMBLE")
```

## Install development version from GitHub

```r
devtools::install_github("LucasFett/RUMBLE")
```

---

# Quick Start

```r
library(RUMBLE)
library(phyloseq)

# Load example data
ps <- readRDS(system.file(
    "extdata",
    "PRJEB38465_phyloseq_com_metadados_completos.rds",
    package = "RUMBLE"
))

# Run the full pipeline
set.seed(42)

results <- RUMBLE(
    input = ps,
    outcome_var = "ses",
    class_of_interest = "Low",  # Mandatory: specify your target condition
    tax_level = "Genus",
    output_dir = "my_results"
)

# View results
print(results$plots$roc)
print(results$metrics)
```

---

# Features

- Accepts both **phyloseq objects** and **count matrices** as input  
- Automated preprocessing:
  - taxonomic aggregation  
  - prevalence filtering  
  - compositional normalization (**CLR** or **rCLR**)  
- Four ML algorithms with hyperparameter tuning via parallelized cross-validation  
- Optional `metric_cutoffs` to filter out underperforming models before SHAP interpretation  
- SHAP-based consensus biomarker ranking using Spearman correlation to prevent sparsity-induced directionality errors  
- Publication-ready visualizations:
  - ROC curves  
  - Confusion matrices  
  - SHAP beeswarm plots  
  - Integrated biomarker prevalence  

---

# ⚠️ Experimental Feature: Robust CLR (rCLR) Normalization

By default, RUMBLE uses standard CLR transformation. However, CLR relies on artificial pseudocounts, transforming biological "absence" (true zeros) into fluctuating negative values based on each patient's geometric mean. This can introduce structural noise, especially for tree-based models like **Random Forest** and **XGBoost**.

To address this, RUMBLE offers the **Robust CLR (rCLR)** option:

```r
normalization_method = "rclr"
```

This approach calculates the transformation strictly on non-zero values, preserving the true zero structure of the sparse matrix.

---

# Important Methodological Limitations of rCLR

While rCLR shields tree-based models from pseudocount noise, it is currently considered experimental in RUMBLE due to the following mathematical biases:

### 1. Model Misinterpretation
Linear (**Elastic Net**) and distance-based (**KNN**) models may misinterpret true zeros (absence) as intermediate values between low and high abundance.

### 2. Directionality Disruption
The default SHAP direction metric (Spearman correlation) evaluates continuous ranks and struggles to capture the non-monotonic biological trajectory:

```
Absence → Low Abundance → High Abundance
```

### 3. SHAP Deflation
Without pseudocount inflation, true zeros (0) reduce the global SHAP impact of sparse taxa, potentially removing biologically relevant biomarkers from top rankings.

---

We are actively researching solutions to improve SHAP directionality under rCLR for future releases.

---

# Citation

If you use RUMBLE in your research, please cite:

Fett, L. (2026). *RUMBLE: Robust Unified Microbiome Biomarker Learning Engine*. R package version 0.99.0.

---

# License

GPL (>= 3)
