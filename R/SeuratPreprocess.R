#' Create a Seurat object from count data
#'
#' This function creates a Seurat object from count data using Seurat::CreateSeuratObject.
#'
#' @export
#' @importFrom Seurat CreateSeuratObject
#' @param counts_data A matrix or data frame of count data.
#' @param ... Additional arguments to be passed to Seurat::CreateSeuratObject.
#' @return A Seurat object.
SeuratPreprocess <- function(counts_data, ...) {
  seurat_object <- Seurat::CreateSeuratObject(counts = counts_data, project = "project_title", min.cells = 3, min.features = 200)
}
