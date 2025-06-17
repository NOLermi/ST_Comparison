# Xenium Data analysis for unimodal and multi-modal versions

library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)

# Create Seurat object
path <- "/path to data/"
MESO1<- LoadXenium(path, fov = "fov",assay = "Xenium")

# Remove cells with fewer than 10 transcript counts
MESO1<- subset(x=MESO1,subset = nCount_Xenium > 10))

## Add slide ID
length(colnames(MESO1))
Slide_ID <-rep(c("MESO1"),times=length(colnames(MESO1))
Slide_ID<-data.frame(Slide_ID)
Slide_ID$idents <- Idents(MESO1)
MESO1<- AddMetaData(MESO1, metadata=Slide_ID,col.name = "Slide_ID")


## MESO2 

path <- "/path to data/"
MESO2<- LoadXenium(path, fov = "fov",assay = "Xenium")

# Remove cells with fewer than 10 transcript counts
MESO2<- subset(x=MESO2,subset = nCount_Xenium > 10))

## Add slide ID
length(colnames(MESO2))
Slide_ID <-rep(c("MESO2"),times=length(colnames(MESO2))
Slide_ID<-data.frame(Slide_ID)
Slide_ID$idents <- Idents(MESO2)
MESO2<- AddMetaData(MESO1, metadata=Slide_ID,col.name = "Slide_ID")

               
# MErge raw data 
MESO<- merge(MESO1, y = MESO2, add.cell.ids = c("MESO1", "MESO2"), project = "MESO")
MESO
head(slot(object = MESO, name = "meta.data"))

MESO<- Seurat::NormalizeData(MESO,normalization.method = "LogNormalize",verbose = T,assay = "Xenium")
MESO <- FindVariableFeatures(MESO,selection.method = "vst", nfeatures = 1000) 
MESO<- ScaleData(MESO,verbose = T)

MESO <- RunPCA(MESO, npcs = 30, features = rownames(MESO))

## Integrate with HARMONY
library(harmony)
MESO<-MESO %>% 
  RunHarmony("Slide_ID", plot_convergence = TRUE)
MESO <- MESO %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()


# Find cell markers per clusters
MESO <- JoinLayers(MESO)
MESO_markers <- FindAllMarkers(MESO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

saveRDS(MESO,file="MESO_Xenium.rds")
