# CosMx Data analysis

library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)

# Create Seurat object
Icon1<- LoadNanostring(data.dir = "/path to Cosmx data/",
                          fov="Icon1", assay = "Nanostring")

# Remove cells with fewer than 30 transcript counts
Icon1<- subset(x=Icon1,subset = nCount_Nanostring > 30)

# (Remove cells with areas more than five times the geometric mean area of all cells. 
# Geometric mean area of all cells.

Geomean <- exp(mean(log(Icon1$Area)))
Icon1<-  subset(Icon1, subset = Area < 5*Geomean)

# Check length(colnames(Icon1)) to find times value
Slide_ID <- rep(c("Icon1"),times=length(colnames(Icon1)))
names(Slide_ID) <- colnames(x = Icon1)

# Add Slide_ID column to the metadata
Icon1 <- AddMetaData(Icon1,Slide_ID, col.name = "slide_ID" )

# Create Seurat object using background subtracted matrix
Icon2<- LoadNanostring(data.dir = "/path to Cosmx data/",
                          fov="Icon2", assay = "Nanostring")

# Remove cells with fewer than 30 transcript counts
Icon2<- subset(x=Icon1,subset = nCount_Nanostring > 30)

# (Remove cells with areas more than five times the geometric mean area of all cells. 
# Geometric mean area of all cells.

Geomean <- exp(mean(log(Icon2$Area)))
Icon2<-  subset(Icon2, subset = Area < 5*Geomean)

# Check length(colnames(Icon2)) to find times value
Slide_ID <- rep(c("Icon2"),times=length(colnames(Icon2)))
names(Slide_ID) <- colnames(x = Icon2)

# Add Slide_ID column to the metadata
Icon2 <- AddMetaData(Icon2,Slide_ID, col.name = "slide_ID" )

# Merge datasets
library(harmony)

ICON<- merge(Icon1, y = Icon2, add.cell.ids = c("Icon1", "Icon2"), project = "ICON")
ICON
head(slot(object = ICON, name = "meta.data"))

unique(ICON$slide_ID)
# create a list of the original data that we loaded to start with
ICON <- Seurat::NormalizeData(ICON,normalization.method = "LogNormalize",verbose = T)
ICON <- FindVariableFeatures(ICON,selection.method = "vst", nfeatures = 1000) 
ICON<- ScaleData(ICON,verbose = T)
ICON <- RunPCA(ICON,npcs = 30, verbose = T)


## Integrate with HARMONY
ICON <- ICON %>% 
  RunHarmony("slide_ID", plot_convergence = TRUE)

ICON <- ICON %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

# Find cell markers per clusters
ICON <- JoinLayers(ICON)
ICON_markers <- FindAllMarkers(ICON, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

saveRDS(ICON,file="ICON_CosMx.rds")
