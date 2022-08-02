library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(gridExtra)

setwd('~/Box/RWorkSpace/Spatial_final/Spatial/Rfiles/')

# Read in data (specify ".h5" data file path)
data = Load10X_Spatial(data.dir = "../A2_KPx2_6min/outs/")
data$orig.ident <- 'A2'

outliers <- read.csv('../A2_KPx2_6min/outs/A2_Outliers.csv')
rownames(outliers) <- outliers$Barcode
outliers <- outliers[colnames(data),]
data$outliers <- outliers$Outliers
data <- subset(data, outliers == 'Outliers', invert = T)

plot1 <- VlnPlot(data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(data, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

data <- SCTransform(data, assay = "Spatial", method = "glmGamPoi", verbose = FALSE)

SpatialFeaturePlot(data, features = c("Sftpc", "Scgb1a1"))

data <- RunPCA(data, assay = "SCT", verbose = FALSE)
data <- RunUMAP(data, reduction = "pca", dims = 1:50)
data <- FindNeighbors(data, reduction = "pca", dims = 1:50)
data <- FindClusters(data, resolution = 0.5)

p1 <- DimPlot(data, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(data, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(data, cells.highlight = CellsByIdentities(object = data, idents = c(0:7)), facet.highlight = TRUE, ncol = 4)

load('../../Integration/T_cells/transferred_rna_reference.RData')

anchors <- FindTransferAnchors(reference = rna_reference, query = data, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = rna_reference$cell_type, prediction.assay = TRUE, 
                                  weight.reduction = data[["pca"]], dims = 1:50)
data[["predictions"]] <- predictions.assay

DefaultAssay(data) <- "predictions"

# pdf('spatial_predicted_cell_types.pdf', width = 32, height = 24)
# SpatialFeaturePlot(data, features = c('Alveolar epithelial cells', 'Club cells', 'Fibroblast', 'Endothelial cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'T cells', 'NK cells', 'Erythrocytes'), pt.size.factor = 1.6, ncol = 4, crop = TRUE)
# dev.off()

DefaultAssay(data) <- "SCT"
data <- FindSpatiallyVariableFeatures(data, assay = "SCT", features = VariableFeatures(data), 
                                      selection.method = "markvariogram")
top.features <- head(SpatiallyVariableFeatures(data, selection.method = "markvariogram"), 6)
SpatialFeaturePlot(data, features = top.features, ncol = 3, alpha = c(0.1, 1))

saveRDS(data, file = 'A2_KPx2_6min.rds')
