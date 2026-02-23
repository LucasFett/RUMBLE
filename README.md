# RUMBLE

**R**obust **U**nified **M**icrobiome **B**iomarker **L**earning **E**ngine

RUMBLE is an R/Bioconductor package that implements a consensus machine
learning framework for identifying microbial biomarkers from
high-throughput sequencing data.

## Overview

RUMBLE integrates multiple classification algorithms (Random Forest,
XGBoost, Elastic Net, and KNN) via the `tidymodels` ecosystem and
provides model-agnostic interpretability through SHAP values and
permutation importance (via `DALEX`). The consensus approach ranks
biomarkers across all models, providing robust and reproducible feature
selection for microbiome research.

## Installation

```r
# Install from Bioconductor (when available)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RUMBLE")

# Install development version from GitHub
devtools::install_github("LucasFett/RUMBLE")
```

## Quick Start

```r
library(RUMBLE)
library(phyloseq)

# Load example data
ps <- readRDS(system.file("extdata",
    "PRJEB38465_phyloseq_com_metadados_completos.rds",
    package = "RUMBLE"))

# Run the full pipeline
set.seed(42)
results <- RUMBLE(
    input = ps,
    outcome_var = "ses",
    tax_level = "Genus",
    output_dir = "my_results"
)

# View results
print(results$plots$roc)
print(results$metrics)
```

## Features

- Accepts both `phyloseq` objects and count matrices as input
- Automated preprocessing: taxonomic aggregation, prevalence filtering,
  CLR transformation
- Four ML algorithms with hyperparameter tuning via cross-validation
- SHAP-based consensus biomarker ranking
- Seven publication-ready visualizations
- Optional removal of unclassified taxa

## Citation

If you use RUMBLE in your research, please cite:

> Fett, L. (2026). RUMBLE: Robust Unified Microbiome Biomarker Learning
> Engine. R package version 0.99.0.

## License

GPL (>= 3)
