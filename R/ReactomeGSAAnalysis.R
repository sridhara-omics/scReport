

#' Title
#'
#' @param lowdim_seurat_object Seurat object that has clusters information
#'
#' @return list
#' @export
#'
#' @importFrom ReactomeGSA analyse_sc_clusters
#' @importFrom ReactomeGSA pathways
ReactomeData <- function(lowdim_seurat_object){
  gsva_result <- NULL
  pathway_expression <- NULL

  gsva_result <- ReactomeGSA::analyse_sc_clusters(lowdim_seurat_object)

  pathway_expression <- ReactomeGSA::pathways(gsva_result)

  colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))

  max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(values), max = max(values)))
  }))

  max_difference$diff <- max_difference$max - max_difference$min
  max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

  return(list(gsva_result, pathway_expression, max_difference))
  #plot_gsva_pathway(gsva_result, pathway_id = rownames(max_difference)[1])
  #plot_gsva_heatmap(gsva_result, max_pathways = 15, margins = c(6,20))

}

