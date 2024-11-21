
## Seurat wrapper to report expressed markers and associated Reactome pathways.

This repository is a simple R package that has 4 main functionalities:  
1. SeuratPreprocess function to transform counts data to scaled seurat object,  
2. SeuratLowDim function to convert the scaled object to low-dimensinoal object with cluster information,  
3. SeuratMarkers function identifies the entire list of markers, along with significant markers (based on minimum percent of cells) and  
4. ReactomeData function to identify the Reactome GSA pathways on the expressed genes in the clusters.  

```{r cars}
library(Seurat)
library(ReactomeGSA)
library(scReport)
```

```{r}
# Loading sample counts data (Mouse cell atlas)
mca.matrix <- readRDS(file = "data/MCA_merged_mat.rds")
mca.metadata <- read.csv(file = "data/MCA_All-batch-removed-assignments.csv", row.names = 1)
mca.matrix.1K <- mca.matrix[,1:1000]
```

[1]. SeuratPreprocess Function to convert Counts Data to Normalized and Scaled Seurat object for highly variable genes.  
```{r SeuratPreprocess function}
scaled_seurat_object <- SeuratPreprocess(mca.matrix.1K)
```

[2]. SeuratLowDim Function to convert Scaled Seurat Object from [1] to object that has clusters identified and data transformed to visualize in 2d (ie, PCA followed by t-SNE and UMAP).  
```{r SeuratLowDim function}
low_dim_object <- SeuratLowDim(scaled_seurat_object)
```

[3]. SeuratMarkers Function to convert seurat object to all markers list, along with significant markers (based on minimum percent of cells in the cluster).  
```{r SeuratMarkers function}
Markers_list <- SeuratMarkers(low_dim_object)
```

[4]. ReactomeData Function to convert seurat object to get the pathways identified using ReactomeGSA R package.  

```{r ReactomeData function}
Reactome_pathways_object <- ReactomeData(low_dim_object)
```
