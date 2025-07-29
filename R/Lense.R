#' Lense: Optimizing Data Preprocessing in Single-Cell Omics Using LLMs
#'
#' Applies multiple preprocessing pipelines to a single-cell Seurat object,
#' generates UMAP plots, and selects the optimal pipeline using GPT-4o-based visual comparison.
#' All UMAP plots are saved automatically. The function returns the best Seurat object
#' corresponding to the optimal preprocessing pipeline.
#'
#' @param seurat_obj A Seurat object after cell filtering and quality control.
#' @param normalization_methods Character vector of normalization methods to evaluate. Default includes
#' \code{c("SCT", "None", "LogNormalize", "CLR", "RC", "LogCount")}.
#' @param nfeatures_vals Character vector specifying feature selection strategies to evaluate.
#' Accepts \code{"HVG"} (select highly variable genes) and \code{"all"} (use all genes).
#' @param pc_vals A numeric vector specifying the number of principal components to evaluate (default: \code{c(5, 20)}).
#' @param resolution_vals Numeric vector of clustering resolution values to try (default: \code{c(0.2, 0.5, 1)}).
#' @param api_key Your OpenAI API key. Defaults to the environment variable \code{OPENAI_API_KEY}.
#' @param output_dir Directory to save UMAP plots. If \code{NULL}, uses a system temporary directory.
#' @param verbose Logical. If \code{TRUE}, print intermediate progress messages. Default is \code{FALSE}.
#'
#' @return A Seurat object corresponding to the best-performing preprocessing pipeline selected by GPT-4o.
#'
#' @import ggplot2
#'
#' @examples
#' \dontrun{
#' Sys.setenv(OPENAI_API_KEY = "your_api_key")
#' best_obj <- Lense(seurat_obj)
#' }
#'
#' @export

Lense <- function(seurat_obj,
                  normalization_methods = c("SCT", "None", "LogNormalize", "LogCount", "RC", "CLR"),
                  nfeatures_vals = c("HVG", "all"),
                  pc_vals = c(5, 20),
                  resolution_vals = c(0.2, 0.5, 1),
                  api_key = Sys.getenv("OPENAI_API_KEY"),
                  output_dir = NULL,
                  verbose = FALSE) {

  if (is.null(output_dir)) {
    output_dir <- file.path(tempdir(), "lense_umaps")
  }
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (verbose) print(paste("All UMAPs will be saved in:", output_dir))

  if (is.null(api_key) || api_key == "" || startsWith(api_key, "your_openai_API_key")) {
    stop("OpenAI API key not found or invalid.")
  }

  seurat_obj <- subset(seurat_obj, subset = is.finite(nCount_RNA) & nCount_RNA > 0)
  if (ncol(seurat_obj) == 0) stop("All cells removed during filtering.")

  I <- nrow(seurat_obj)
  computed_nfeatures <- 2 * 10^(ceiling(log10(I / 2)) - 1)

  DefaultAssay(seurat_obj) <- "RNA"
  pipelines <- list()
  objects_list <- list()

  for (norm_method in normalization_methods) {
    for (nfeatures in nfeatures_vals) {
      for (pc in pc_vals) {
        for (res in resolution_vals) {
          label <- paste0("Norm=", norm_method, " | nFeatures=", nfeatures, " | PCs=", pc, " | Res=", res)
          if (verbose) print(paste("Generating:", label))

          if (norm_method == "SCT") {
            normalized_obj <- SCTransform(seurat_obj, verbose = FALSE)
          } else if (norm_method == "None") {
            normalized_obj <- seurat_obj
            normalized_obj <- ScaleData(normalized_obj, features = rownames(normalized_obj))
          } else if (norm_method == "LogCount") {
            normalized_obj <- seurat_obj
            DefaultAssay(normalized_obj) <- "RNA"
            raw_counts <- GetAssayData(normalized_obj, layer = "counts")
            log_counts <- log1p(raw_counts)
            log_assay <- CreateAssayObject(counts = raw_counts)
            log_assay <- SetAssayData(log_assay, layer = "data", new.data = log_counts)
            normalized_obj[["LogAssay"]] <- log_assay
            DefaultAssay(normalized_obj) <- "LogAssay"
            normalized_obj <- ScaleData(normalized_obj, features = rownames(normalized_obj))
          } else {
            normalized_obj <- NormalizeData(seurat_obj, normalization.method = norm_method)
            normalized_obj <- ScaleData(normalized_obj, features = rownames(normalized_obj))
          }

          if (nfeatures == "all") {
            variable_features <- rownames(normalized_obj)
          } else {
            normalized_obj <- FindVariableFeatures(normalized_obj, nfeatures = computed_nfeatures)
            variable_features <- VariableFeatures(normalized_obj)
          }

          pca_obj <- RunPCA(normalized_obj, npcs = pc, features = variable_features)
          umap_obj <- RunUMAP(pca_obj, dims = 1:pc)
          umap_obj <- FindNeighbors(umap_obj, reduction = "pca", dims = 1:pc)
          umap_obj <- FindClusters(umap_obj, resolution = res)

          umap_plot <- DimPlot(umap_obj, reduction = "umap", group.by = "seurat_clusters")
          plot_filename <- file.path(output_dir, paste0(
            "UMAP_",
            gsub("[=| |\\|]", "_", label),
            ".png"
          ))
          ggsave(filename = plot_filename, plot = umap_plot, width = 10, height = 8)

          pipelines[[length(pipelines) + 1]] <- list(
            name = label,
            file = plot_filename,
            seurat = umap_obj
          )
        }
      }
    }
  }

  current_winner <- pipelines[[1]]
  log_comparisons <- list()

  for (i in 2:length(pipelines)) {
    challenger <- pipelines[[i]]
    if (verbose) print(paste("Comparing:", current_winner$name, "vs", challenger$name))

    if (!file.exists(current_winner$file) || !file.exists(challenger$file)) {
      stop("Image missing before comparison: ", current_winner$file, " or ", challenger$file)
    }

    comparison_result <- gpt_compare_umap(current_winner$file, challenger$file, api_key, verbose = verbose)
    winner <- if (comparison_result$result == 1) current_winner else challenger

    log_comparisons[[length(log_comparisons) + 1]] <- list(
      round = i - 1,
      winner = winner$name,
      loser = if (comparison_result$result == 1) challenger$name else current_winner$name,
      gpt_reply = comparison_result$raw_reply
    )

    current_winner <- winner
  }

  if (verbose) message("Final Selected Pipeline:", current_winner$name)

  return(current_winner$seurat)
}
