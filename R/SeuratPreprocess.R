#' Create a Seurat object from count data
#'
#' This function creates a Seurat object from count data using Seurat::CreateSeuratObject.
#'
#' @export
#' @importFrom Seurat CreateSeuratObject
#' @importFrom Seurat NormalizeData
#' @importFrom Seurat FindVariableFeatures
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom Seurat ScaleData
#' @param counts_data A matrix or data frame of count data.
#' @param ... Additional arguments to be passed to Seurat::CreateSeuratObject.
#' @return A Seurat object.
SeuratPreprocess <- function(counts_data, ...) {
  n_size = dim(mca.matrix.10K)[2]

  seurat_object <- Seurat::CreateSeuratObject(counts = counts_data, project = "project_title", min.cells = 3, min.features = 200)

  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = n_size)

  seurat_object <- FindVariableFeatures(seurat_object)

  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")

  seurat_object <- ScaleData(seurat_object, vars.to.regress = "percent.mt")
}
