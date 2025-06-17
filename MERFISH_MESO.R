

vizgen.input <- ReadVizgen(data.dir = "/path to data/", type = "centroids")
## remove the Blank-* control probes
blank_index<- which(stringr::str_detect(rownames(vizgen.input$transcripts), "^Blank"))

transcripts<-vizgen.input$transcripts[-blank_index, ]

dim(vizgen.input$transcripts)

# Create seurat object
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


vizgen.obj<- subset(x=vizgen.obj,subset = nCount_Vizgen > 10)
