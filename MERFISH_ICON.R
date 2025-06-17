# MERFISH Data analysis 

library(dplyr)
library(ggplot2)
library(patchwork)
library(Seurat)

# Create Seurat object

## remove the Blank-* control probes
blank_index<- which(stringr::str_detect(rownames(vizgen.input$transcripts), "^Blank"))
# blanks are empty
transcripts<-vizgen.input$transcripts

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


