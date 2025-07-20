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


# 🚀 Run the Lense pipeline:
# This will:
# - Apply 24 preprocessing pipelines
# - Generate UMAP plots under each combination
# - Compare them pairwise using GPT-4o (via OpenAI API)
# - Automatically select the optimal preprocessing pipeline

result <- Lense(seurat_obj,
                nfeatures_vals = NULL,                   # NULL = use dynamic formula + all genes
                pc_vals = c(5, 20),                      # Principal Components: typically 5 and 20
                resolution = 0.3,                        # Clustering resolution for Seurat
                output_dir = "/path/to/output/folder"    # Folder to save all generated UMAP plots
)

# ✅ View results:
result$best_pipeline          # Selected optimal preprocessing pipeline
View(result$comparisons_log)  # Full comparison log (optional)
print(result$umap_folder)     # Path to folder with saved UMAP plots

```

### ⚠️ Warning: Do NOT share your API key publicly or upload it to GitHub.

## Introduction
Preprocessing is critical for single-cell omics analyses, but default pipelines may underperform on diverse datasets, especially
from emerging platforms like spatial transcriptomics. We introduce Lense, a language-model-guided method that automatically
selects optimal preprocessing by comparing UMAP plots across pipeline variants. Integrated with Seurat, Lense streamlines
analysis and improves preprocessing robustness without requiring manual tuning.

## Contact

Authors: Jingyun Liu (jingyun.liu2@duke.edu), Zhicheng Ji (zhicheng.ji@duke.edu).

Report bugs and provide suggestions by sending email to the maintainer Dr. Wenpin Hou (jingyun.liu2@duke.edu) or open a new issue on this Github page.


