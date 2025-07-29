# Lense: Optimizing data preprocessing in single-cell omics using LLMs


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
best_obj <- Lense(seurat_obj)

# ✅ Check the result: it's a standard Seurat object
print(best_obj)

# 📊 Visualize the final selected UMAP
DimPlot(best_obj, reduction = "umap", group.by = "seurat_clusters")


```

### ⚠️ Warning: Do NOT share your API key publicly or upload it to GitHub.

## Introduction
Preprocessing is critical for single-cell omics analyses, but default pipelines may underperform on diverse datasets, especially
from emerging platforms like spatial transcriptomics. We introduce Lense, a language-model-guided method that automatically
selects optimal preprocessing by comparing UMAP plots across pipeline variants. Integrated with Seurat, Lense streamlines
analysis and improves preprocessing robustness without requiring manual tuning.

## Contact

Authors: Jingyun Liu (jingyun.liu2@duke.edu), Zhicheng Ji (zhicheng.ji@duke.edu).

Report bugs and provide suggestions by sending email to the maintainer Jingyun Liu (jingyun.liu2@duke.edu) or open a new issue on this Github page.


