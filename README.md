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

remotes::install_github("jingyun20/lense")
```

##  🚀 Quick start

```{r eval = FALSE}
# IMPORTANT! Assign your OpenAI API key
Sys.setenv(OPENAI_API_KEY = 'your_openai_API_key')

# Load packages
library(lense)
library(Seurat)

# Run preprocessing pipelines (Seurat object must be pre-filtered)
# This step generates UMAP plots under 24 combinations of preprocessing parameters
preprocess_pipelines(seurat_obj,
                     nfeatures_vals = c(200, 348),
                     pc_vals = c(5, 20))

# After generating UMAP plots, compare any two using GPT-4o:
res <- compare_umap("path/to/umap1.png", "path/to/umap2.png")

# GPT-4o selects the better UMAP
print(res$result)     # 1 or 2
print(res$raw_reply)  # Full GPT reply

```

### ⚠️ Warning: Do NOT share your API key publicly or upload it to GitHub.

## Introduction
Preprocessing is critical for single-cell omics analyses, but default pipelines may underperform on diverse datasets, especially
from emerging platforms like spatial transcriptomics. We introduce Lense, a language-model-guided method that automatically
selects optimal preprocessing by comparing UMAP plots across pipeline variants. Integrated with Seurat, Lense streamlines
analysis and improves preprocessing robustness without requiring manual tuning.

## Contact

Authors: Jingyun Liu (jingyun.liu2@duke.edu), Zhicheng Ji (zhicheng.ji@duke.edu).


