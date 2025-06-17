# MERFISH Data analysis 

library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)

# Create Seurat object

vizgen.input <- ReadVizgen(data.dir = "/path to data/", type = "centroids")
## remove the Blank-* control probes
blank_index<- which(stringr::str_detect(rownames(vizgen.input$transcripts), "^Blank"))

transcripts<-vizgen.input$transcripts[-blank_index, ]

dim(vizgen.input$transcripts)


vizgen.obj<- CreateSeuratObject(counts = transcripts, assay = "Vizgen")

cents <- CreateCentroids(vizgen.input$centroids)
segmentations.data <- list(
  "centroids" = cents,
  "segmentation" = NULL
)

coords <- CreateFOV(
  coords = segmentations.data,
  type = c("segmentation", "centroids"),
  molecules = NULL,
  assay = "Vizgen"
)

vizgen.obj[["vizgen.obj"]] <- coords

# Remove cells with fewer than 10 transcripts
vizgen.obj<- subset(x=vizgen.obj,subset = nCount_Vizgen > 10)

vizgen.obj <- NormalizeData(vizgen.obj, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  ScaleData() 
vizgen.obj <- RunPCA(vizgen.obj, npcs = 30, features = rownames(vizgen.obj))
vizgen.obj <- RunUMAP(vizgen.obj, dims = 1:30)
vizgen.obj <- FindNeighbors(vizgen.obj, reduction = "pca", dims = 1:30)
vizgen.obj <- FindClusters(vizgen.obj, resolution = 0.8)

saveRDS(vizgen.obj, file ="ICON_MERFISH.rds")
