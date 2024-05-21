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
#' @importFrom Seurat RunPCA
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat RunTSNE
#' @importFrom Seurat RunUMAP
#' @param counts_data A matrix or data frame of count data.
#' @param ... Additional arguments to be passed to Seurat::CreateSeuratObject.
#' @return A Seurat object.
SeuratPreprocess <- function(counts_data, ...) {
  n_size <- dim(counts_data)[2]

  n_pcs <- round(n_size/100)

  n_dims <- round(n_pcs/1.25)

  seurat_object <- Seurat::CreateSeuratObject(counts = counts_data, project = "project_title", min.cells = 3, min.features = 200)

  seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = n_size)

  seurat_object <- FindVariableFeatures(seurat_object)

  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")

  seurat_object <- ScaleData(seurat_object, vars.to.regress = "percent.mt")

  seurat_object <- RunPCA(seurat_object, npcs = n_pcs, ndims.print = 1:5, nfeatures.print = 5)

  seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:n_dims, nn.eps = 0.5)

  seurat_object <- FindClusters(seurat_object, resolution = 3, n.start = 10)

  seurat_object <- RunTSNE(seurat_object, dims = 1:n_dims)

  seurat_object <- RunUMAP(seurat_object, dims = 1:n_dims, min.dist = 0.75)
}
