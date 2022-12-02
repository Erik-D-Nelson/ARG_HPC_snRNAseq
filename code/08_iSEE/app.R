library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")

load("sce_iSEE.rda", verbose = TRUE)

source("initial.R", print.eval = TRUE)

colData(sce) <- cbind(
  colData(sce)[, !colnames(colData(sce)) %in% c("Sample", "cellType")],
  colData(sce)[, c("cellType", "Sample")]
)

sce$Sample <- as.factor(sce$Sample)

sce <- registerAppOptions(sce, color.maxlevels = 18)
iSEE(
    sce,
    appTitle = "2022_HPC_ARG",
    initial = initial,
    colormap = ExperimentColorMap(colData = list(
        Sample = function(n) {
            cols <- paletteer::paletteer_d(
                palette = "RColorBrewer::Dark2",
                n = length(unique(sce$Sample))
            )
            cols <- as.vector(cols)
            names(cols) <- levels(sce$Sample)
            return(cols)
        },
         cellType = function(n) {
             cols <- paletteer::paletteer_d(
                 palette = "Polychrome::palette36",
                 n = length(levels(sce$cellType))
             )
             cols <- as.vector(cols)
             names(cols) <- levels(sce$cellType)
             return(cols)
         }
    ))
)
