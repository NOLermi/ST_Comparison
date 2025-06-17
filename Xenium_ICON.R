# Xenium Data analysis for unimodal and multi-modal versions

library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)

# Create Seurat object
path <- "/path to data/"
ICON<- LoadXenium(path, fov = "fov",assay = "Xenium")

# Remove cells with fewer than 10 transcript counts
ICON<- subset(x=ICON,subset = nCount_Xenium > 10))

ICON <- NormalizeData(ICON, assay="Xenium")
ICON <-ScaleData(ICON, assay="Xenium")
ICON <- RunPCA(ICON, npcs = 30, features = rownames(ICON))
ICON <- RunUMAP(ICON dims = 1:30)
ICON <- FindNeighbors(ICON, reduction = "pca", dims = 1:30)
ICON <- FindClusters(ICON, resolution = 0.8,future.seed = NULL)

# Find cell markers per clusters
ICON <- JoinLayers(ICON)
ICON_markers <- FindAllMarkers(ICON, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

saveRDS(ICON,file="ICON_Xenium.rds")
