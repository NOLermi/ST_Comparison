# CosMx Data analysis

library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)

# Create Seurat object
MESO1<- LoadNanostring(data.dir = "/path to Cosmx data/",
                          fov="MESO1", assay = "Nanostring")

# Remove cells with fewer than 30 transcript counts
MESO1<- subset(x=MESO1,subset = nCount_Nanostring > 30)

# (Remove cells with areas more than five times the geometric mean area of all cells. 
# Geometric mean area of all cells.

Geomean <- exp(mean(log(MESO1$Area)))
MESO1<-  subset(MESO1, subset = Area < 5*Geomean)

# Check length(colnames(MESO1)) to find times value
Slide_ID <- rep(c("MESO1"),times=length(colnames(MESO1)))
names(Slide_ID) <- colnames(x = MESO1)

# Add Slide_ID column to the metadata
MESO1 <- AddMetaData(MESO1,Slide_ID, col.name = "slide_ID" )

# Create Seurat object using background subtracted matrix
MESO2<- LoadNanostring(data.dir = "/path to Cosmx data/",
                          fov="MESO2", assay = "Nanostring")

# Remove cells with fewer than 30 transcript counts
MESO2<- subset(x=MESO1,subset = nCount_Nanostring > 30)

# (Remove cells with areas more than five times the geometric mean area of all cells. 
# Geometric mean area of all cells.

Geomean <- exp(mean(log(MESO2$Area)))
MESO2<-  subset(MESO2, subset = Area < 5*Geomean)

# Check length(colnames(MESO2)) to find times value
Slide_ID <- rep(c("MESO2"),times=length(colnames(MESO2)))
names(Slide_ID) <- colnames(x = MESO2)

# Add Slide_ID column to the metadata
MESO2 <- AddMetaData(MESO2,Slide_ID, col.name = "slide_ID" )

# Merge datasets
library(harmony)

MESO<- merge(MESO1, y = MESO2, add.cell.ids = c("MESO1", "MESO2"), project = "MESO")
MESO
head(slot(object = MESO, name = "meta.data"))

unique(MESO$slide_ID)
# create a list of the original data that we loaded to start with
MESO <- Seurat::NormalizeData(MESO,normalization.method = "LogNormalize",verbose = T)
MESO <- FindVariableFeatures(MESO,selection.method = "vst", nfeatures = 1000) 
MESO<- ScaleData(MESO,verbose = T)
MESO <- RunPCA(MESO,npcs = 30, verbose = T)


## Integrate with HARMONY
MESO <- MESO %>% 
  RunHarmony("slide_ID", plot_convergence = TRUE)

MESO <- MESO %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

# Find cell markers per clusters
MESO <- JoinLayers(MESO)
MESO_markers <- FindAllMarkers(MESO, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

saveRDS(MESO,file="MESO_CosMx.rds")
