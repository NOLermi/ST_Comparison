# CosMx Data analysis

library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)

# Create Seurat object using background subtracted matrix
Icon1<- LoadNanostring(data.dir = "/path to Cosmx data/",
                          fov="Icon1", assay = "Nanostring")

# Remove cells with fewer than 30 transcript counts
Icon1<- subset(x=Icon1,subset = nCount_Nanostring > 30)

# (Remove cells with areas more than five times the geometric mean area of all cells. 
# Geometric mean area of all cells.

Geomean <- exp(mean(log(Icon1$Area)))
Icon1<-  subset(Icon1, subset = Area < 5*Geomean)

# Create Seurat object using background subtracted matrix
Icon2<- LoadNanostring(data.dir = "/path to Cosmx data/",
                          fov="Icon2", assay = "Nanostring")

# Remove cells with fewer than 30 transcript counts
Icon2<- subset(x=Icon1,subset = nCount_Nanostring > 30)

# (Remove cells with areas more than five times the geometric mean area of all cells. 
# Geometric mean area of all cells.

Geomean <- exp(mean(log(Icon2$Area)))
Icon2<-  subset(Icon2, subset = Area < 5*Geomean)
