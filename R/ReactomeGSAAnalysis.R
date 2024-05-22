

#' Title
#'
#' @param lowdim_seurat_object
#'
#' @return list
#' @export
#'
#' @importFrom ReactomeGSA analyse_sc_clusters
#' @importFrom ReactomeGSA pathways
ReactomeData <- function(lowdim_seurat_object){
  gsva_result <- ReactomeGSA::analyse_sc_clusters(lowdim_seurat_object, verbose = TRUE)
  pathway_expression <- ReactomeGSA::pathways(gsva_result)
  return(list(gsva_result, pathway_expression))
  #plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1])
  #plot_gsva_heatmap(gsva_result, max_pathways = 15, margins = c(6,20))

}

