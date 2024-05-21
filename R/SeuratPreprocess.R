library(Seurat)

SeuratPreprocess <- function(counts_data) {
  seurat_object <- Seurat::CreateSeuratObject(counts = counts_data, project = "project_title", min.cells = 3, min.features = 200)
  return(seurat_object)
}
