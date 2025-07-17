#' Generate UMAP Plots Under 24 Preprocessing Pipelines
#'
#' This function applies 24 combinations of preprocessing pipelines to a Seurat object,
#' generates UMAP plots, and saves each as a PNG image.
#'
#' @param seurat_obj Seurat object: input single-cell omics data after cell filtering.
#' @param nfeatures_vals Vector of numbers of highly variable features to select.
#' @param pc_vals Vector of numbers of principal components to use. Default: c(5, 20).
#' @param resolution Resolution parameter for clustering. Default: 0.3.
#'
#' @return Saves UMAP PNG files to the working directory. Returns invisibly.
#' @export
#'
#' @import Seurat ggplot2 Matrix
#'
#' @examples
#' \dontrun{
#' preprocess_pipelines(seurat_obj, c(200, 367), c(5, 20))
#' }
preprocess_pipelines <- function(seurat_obj, nfeatures_vals, pc_vals = c(5, 20), resolution = 0.3) {

  normalization_methods <- list("SCT", "None", "LogNormalize", "CLR", "RC", "LogCount")

  DefaultAssay(seurat_obj) <- "RNA"

  for (norm_method in normalization_methods) {
    for (nfeatures in nfeatures_vals) {
      for (pc in pc_vals) {

        message(sprintf("Processing: Norm=%s | Features=%d | PCs=%d", norm_method, nfeatures, pc))

        # === Normalization Step ===
        if (norm_method == "SCT") {
          normalized_obj <- SCTransform(seurat_obj, verbose = FALSE)

        } else if (norm_method == "None") {
          normalized_obj <- seurat_obj
          normalized_obj <- ScaleData(normalized_obj, features = rownames(normalized_obj))

        } else if (norm_method == "LogCount") {
          normalized_obj <- seurat_obj
          DefaultAssay(normalized_obj) <- "RNA"

          raw_counts <- GetAssayData(normalized_obj, slot = "counts")
          log_counts <- log1p(raw_counts)

          log_assay <- CreateAssayObject(counts = raw_counts, data = log_counts)
          normalized_obj[["LogAssay"]] <- log_assay
          DefaultAssay(normalized_obj) <- "LogAssay"

          normalized_obj <- ScaleData(normalized_obj, features = rownames(normalized_obj))

        } else {
          normalized_obj <- NormalizeData(seurat_obj, normalization.method = norm_method)
          normalized_obj <- ScaleData(normalized_obj, features = rownames(normalized_obj))
        }

        normalized_obj <- FindVariableFeatures(normalized_obj, nfeatures = nfeatures)
        variable_features <- VariableFeatures(normalized_obj)

        pca_obj <- RunPCA(normalized_obj, npcs = pc, features = variable_features)
        umap_obj <- RunUMAP(pca_obj, dims = 1:pc)
        umap_obj <- FindNeighbors(umap_obj, reduction = "pca", dims = 1:pc)
        umap_obj <- FindClusters(umap_obj, resolution = resolution)

        umap_plot <- DimPlot(umap_obj, reduction = "umap", group.by = "seurat_clusters")
        plot_filename <- paste0("UMAP_Norm_", norm_method, "_nFeatures_", nfeatures, "_PCs_", pc, ".png")
        ggsave(filename = plot_filename, plot = umap_plot, width = 10, height = 8)

        message(sprintf("Saved %s", plot_filename))

      }
    }
  }

  invisible(NULL)
}
