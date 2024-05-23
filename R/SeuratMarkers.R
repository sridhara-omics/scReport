#' A thresholded markers list for better calculation of DE genes
#'
#' @param lowdim_seurat_object Seurat object with cluster information
#' @importFrom Seurat FindAllMarkers
#' @return A list
#' @export
#'
#'
SeuratMarkers <- function(lowdim_seurat_object){
  markers_list <- Seurat::FindAllMarkers(lowdim_seurat_object)
  markers_list_threshold <- Seurat::FindAllMarkers(lowdim_seurat_object, min.pct = 0.1)
  return(list(markers_list, markers_list_threshold))
}
