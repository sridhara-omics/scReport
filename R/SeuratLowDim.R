#' Create a Low dimensional Seurat object from scaled seurat object
#'
#' This function converts the transformed data to low-dimensional data for downstream analysis.
#'
#' @export
#' @importFrom Seurat RunPCA
#' @importFrom Seurat FindNeighbors
#' @importFrom Seurat FindClusters
#' @importFrom Seurat RunTSNE
#' @importFrom Seurat RunUMAP
#' @param scaled_seurat_object A scaled Seurat object.
#' @param ... Additional arguments to be passed to convert to low dimensional data frame.
#' @return A Seurat object.
SeuratLowDim <- function(scaled_seurat_object, ...) {
  n_size <- dim(scaled_seurat_object@meta.data)[1]

  n_pcs <- round(n_size/100)

  n_dims <- round(n_pcs/1.25)

  seurat_object <- RunPCA(scaled_seurat_object, npcs = n_pcs, ndims.print = 1:5, nfeatures.print = 5)

  seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:n_dims, nn.eps = 0.5)

  seurat_object <- FindClusters(seurat_object, resolution = 3, n.start = 10)

  seurat_object <- RunTSNE(seurat_object, dims = 1:n_dims)

  seurat_object <- RunUMAP(seurat_object, dims = 1:n_dims, min.dist = 0.75)
}
