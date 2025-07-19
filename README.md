# Lense: LLM-powered preprocessing optimization for single-cell omics


## Installation 

To install the latest version of `Lense` from GitHub, run the following commands in R:

```{r eval = FALSE}
install.packages("Seurat")
install.packages("ggplot2")
install.packages("Matrix")
install.packages("httr")
install.packages("jsonlite")
install.packages("base64enc")

remotes::install_github("jingyun20/Lense")
```

##  🚀 Quick start

```{r eval = FALSE}
# IMPORTANT! Assign your OpenAI API key
Sys.setenv(OPENAI_API_KEY = 'your_openai_API_key')

# Load packages
library(Lense)
library(Seurat)

# ⚙️ Prepare a Seurat object (must be pre-filtered)
# Example: after QC and cell filtering

# 🚀 Run the Lense pipeline:
# This will generate UMAP plots under all preprocessing pipelines,
# compare them pairwise using GPT-4o,
# and automatically select the optimal preprocessing pipeline.

result <- Lense(seurat_obj,
                nfeatures_vals = NULL,  # NULL = use dynamic + all genes
                pc_vals = c(5, 20),     # PCs: typically 5 or 20
                resolution = 0.3）       # Clustering resolutio


```

### ⚠️ Warning: Do NOT share your API key publicly or upload it to GitHub.

## Introduction
Preprocessing is critical for single-cell omics analyses, but default pipelines may underperform on diverse datasets, especially
from emerging platforms like spatial transcriptomics. We introduce Lense, a language-model-guided method that automatically
selects optimal preprocessing by comparing UMAP plots across pipeline variants. Integrated with Seurat, Lense streamlines
analysis and improves preprocessing robustness without requiring manual tuning.

## Contact

Authors: Jingyun Liu (jingyun.liu2@duke.edu), Zhicheng Ji (zhicheng.ji@duke.edu).


