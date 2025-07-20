#' Lense: Automated Selection of Preprocessing Pipeline via GPT-4o or Simulated Mode
#'
#' Applies multiple preprocessing pipelines to a single-cell Seurat object,
#' generates UMAP plots, and selects the optimal pipeline using either GPT-4o-based visual comparison
#' or simulated mode. All UMAP plots are saved automatically. Returns the selected best pipeline
#' and the full pairwise comparison log.
#'
#' @param seurat_obj A Seurat object after cell filtering and quality control.
#' @param nfeatures_vals A vector or list specifying the number of highly variable features (e.g., c(200, 500)) to evaluate.
#' Use \code{NULL} (default) to apply automated selection: using all features and a dynamic value computed as
#' \eqn{a = 2 \times 10^{\lceil \log_{10}(I / 2) \rceil - 1}}, where \eqn{I} is the total number of features (genes).
#' @param pc_vals A numeric vector specifying the number of principal components to evaluate (default: \code{c(5, 20)}).
#' @param resolution Clustering resolution for Seurat clustering (default: \code{0.3}).
#' @param api_key Your OpenAI API key as a character string. Defaults to the environment variable \code{OPENAI_API_KEY}.
#'
#' @return A list containing:
#' \item{best_pipeline}{A character string describing the selected optimal preprocessing pipeline.}
#' \item{comparisons_log}{A list containing the full log of pairwise comparisons.}
#' \item{umap_folder}{A character string with the path to the folder containing all UMAP plots.}
#'
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' \dontrun{
#' Sys.setenv(OPENAI_API_KEY = "your_api_key")
#' result <- Lense(seurat_obj)
#' result$best_pipeline
#' View(result$comparisons_log)
#' }
#'
#'
#' @export

Lense <- function(seurat_obj, nfeatures_vals = NULL, pc_vals = c(5, 20), resolution = 0.3,
                  api_key = Sys.getenv("OPENAI_API_KEY")) {

  output_dir <- file.path(normalizePath(getwd()), "lense_umaps")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  message("All UMAPs will be saved in: ", output_dir)

  if (is.null(api_key) || api_key == "" || startsWith(api_key, "your_openai_API_key")) {
    stop("OpenAI API key not found or invalid.")
  }

  message("Checking and filtering cells with invalid nCount_RNA...")
  seurat_obj <- subset(seurat_obj, subset = is.finite(nCount_RNA) & nCount_RNA > 0)
  if (ncol(seurat_obj) == 0) stop("All cells removed during filtering.")

  message("Using GPT-4o for UMAP visual comparison...")

  I <- nrow(seurat_obj)
  computed_nfeatures <- 2 * 10^(ceiling(log10(I / 2)) - 1)
  message("Total features (I): ", I)
  message("Computed nFeatures based on formula: ", computed_nfeatures)

  if (is.null(nfeatures_vals)) {
    nfeatures_vals <- list("all", computed_nfeatures)
  }

  normalization_methods <- list("SCT", "None", "LogNormalize", "CLR", "RC", "LogCount")
  DefaultAssay(seurat_obj) <- "RNA"
  pipelines <- list()

  for (norm_method in normalization_methods) {
    for (nfeatures in nfeatures_vals) {
      for (pc in pc_vals) {
        message(sprintf("Generating: Norm=%s | Features=%s | PCs=%d",
                        norm_method, as.character(nfeatures), pc))

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
          normalized_obj <- FindVariableFeatures(normalized_obj, nfeatures = nfeatures)
          variable_features <- VariableFeatures(normalized_obj)
        }

        pca_obj <- RunPCA(normalized_obj, npcs = pc, features = variable_features)
        umap_obj <- RunUMAP(pca_obj, dims = 1:pc)
        umap_obj <- FindNeighbors(umap_obj, reduction = "pca", dims = 1:pc)
        umap_obj <- FindClusters(umap_obj, resolution = resolution)

        umap_plot <- DimPlot(umap_obj, reduction = "umap", group.by = "seurat_clusters")
        plot_filename <- file.path(output_dir, paste0(
          "UMAP_",
          gsub("[=| ]", "_", paste0("Norm=", norm_method, "_nFeatures=", nfeatures, "_PCs=", pc)),
          ".png"
        ))

        ggsave(filename = plot_filename, plot = umap_plot, width = 10, height = 8)

        pipelines[[length(pipelines) + 1]] <- list(
          name = paste0("Norm=", norm_method, " | nFeatures=", nfeatures, " | PCs=", pc),
          file = plot_filename
        )
      }
    }
  }

  current_winner <- pipelines[[1]]
  log_comparisons <- list()

  for (i in 2:length(pipelines)) {
    challenger <- pipelines[[i]]
    message(sprintf("Comparing:\n  Winner: %s\n  Challenger: %s", current_winner$name, challenger$name))

    if (!file.exists(current_winner$file) || !file.exists(challenger$file)) {
      stop("Image missing before comparison: ", current_winner$file, " or ", challenger$file)
    }

    comparison_result <- gpt_compare_umap(current_winner$file, challenger$file, api_key)
    winner <- if (comparison_result$result == 1) current_winner else challenger

    log_comparisons[[length(log_comparisons) + 1]] <- list(
      round = i - 1,
      winner = winner$name,
      loser = if (comparison_result$result == 1) challenger$name else current_winner$name,
      gpt_reply = comparison_result$raw_reply
    )

    current_winner <- winner
  }

  message(sprintf("Final Selected Pipeline:\n%s", current_winner$name))

  return(
    best_pipeline = current_winner$name
  )
}
