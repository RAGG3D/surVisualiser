# surVisualiser

Multi-dimensional survival analysis and visualization for TCGA
transcriptome data, with circular-analysis-safe hallmark pathway
scoring for secretome research.

## Overview

**surVisualiser** provides an end-to-end workflow for TCGA cancer
genomics, covering:

| Module | Functions |
|--------|-----------|
| **Expression** | `TCGA_transcript()` -- read & aggregate transcript-level counts |
| **Clinical** | `clinical_combine()` -- merge expression with survival metadata |
| **Immune deconvolution** | `CIBERSORT()` -- cell-type proportion estimation via tidybulk |
| **Survival analysis** | `ggsurvp_cat_item()`, `split_median()`, `strata_p()`, `pw_diff_selection()` |
| **Hallmark GSEA** | `gsea_test_scrto_gene()` -- per-sample pre-ranked GSEA with secretome gene exclusion |
| **Correlation** | `corr_pathway_scrto()` -- Pearson correlation between hallmark scores and secretome expression |
| **Visualization** | `three_heatmaps_of_pth_scrto_TN()` -- Tumor / Normal / T-N difference heatmaps |

### Bundled data

- `secretome_genes` -- curated secreted protein gene list (HPA, VerSeDa)
- `my_ref` -- CIBERSORT reference matrix
- `cellname` -- cell-type name mapping table

## Installation

```r
# Install from Bioconductor (when accepted)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("surVisualiser")

# Or install the development version from GitHub
BiocManager::install("RAGG3D/surVisualiser")
```

Optional packages for the GSEA and heatmap workflow:

```r
BiocManager::install(c(
    "msigdbr", "clusterProfiler",
    "tidyHeatmap", "ComplexHeatmap",
    "circlize", "RColorBrewer"
))
```

## Quick start

### 1. Load expression and clinical data

```r
library(surVisualiser)

brca_tumor <- TCGA_transcript("BRCA", data_dir = "/path/to/data",
                               sample_type = "tumor")

brca_clinical <- clinical_combine(brca_tumor, cancer = "BRCA",
                                   clinical_dir = "/path/to/clinical")
```

### 2. Survival analysis

```r
surv_data <- split_median(brca_clinical, cat = "cancer_stage", item = "TP53")

ggsurvp_cat_item(surv_data, p = TRUE, nrow = 2,
                 palette = c("#E41A1C", "#377EB8"),
                 xlab = "Days", ylab = "Survival Probability",
                 title = "TP53 and Survival", size0 = 12)
```

### 3. Secretome-aware hallmark scoring

```r
data(secretome_genes)

# GSEA with secretome genes excluded from ranking
hallmarks <- gsea_test_scrto_gene(
    data = brca_tumor,
    rank_by = "raw_count_scaled_adjusted",
    excluded_genes = secretome_genes$symbol
)

# Pearson correlation: hallmark scores vs secretome expression
corr_tumor <- corr_pathway_scrto(
    hallmark_scores = hallmarks,
    expression = brca_tumor,
    secretome_symbols = secretome_genes$symbol
)
```

### 4. Three-panel heatmap

```r
# Run step 3 for both tumor and normal, then:
three_heatmaps_of_pth_scrto_TN(
    tumor_corr = corr_tumor,
    normal_corr = corr_normal,
    cancer = "BRCA",
    output_dir = "figures/BRCA"
)
```

## Why exclude secretome genes?

When hallmark gene sets contain secretome genes and the downstream
analysis correlates secretome expression against hallmark scores, the
correlation is inflated because the same genes contribute to both
sides. `gsea_test_scrto_gene()` removes the specified gene list from
the expression matrix before ranking, so the resulting enrichment
scores are independent of the gene set being correlated against.

## License

Artistic-2.0
