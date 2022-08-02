library(RCTD)
library(Matrix)
library(Seurat)
library(ComplexHeatmap)

setwd('~/Box/RWorkSpace/Spatial_combine/Deconvolution/')

read_spatial_RCTD <- function(spatial, rctd){
  data <- readRDS(spatial)
  load(rctd)
  results <- myRCTD@results
  norm_weights = as.matrix(sweep(results$weights, 1, rowSums(results$weights), '/'))
  norm_weights <- t(norm_weights)
  ncol <- length(colnames(data)[!(colnames(data) %in% colnames(norm_weights))])
  if (ncol > 0) {
    tmp <- matrix(0, nrow = nrow(norm_weights), ncol = ncol)
    colnames(tmp) <- colnames(data)[!(colnames(data) %in% colnames(norm_weights))]
    norm_weights <- cbind(norm_weights, tmp)[,colnames(data)]
  }
  data[['RCTD']] <- CreateAssayObject(data = norm_weights)
  DefaultAssay(data) <- "RCTD"
  pdf(gsub('.RData', '.pdf', rctd), width = 16, height = 12)
  print(SpatialFeaturePlot(data, features = rownames(data), pt.size.factor = 1.6, ncol = 4, crop = TRUE))
  dev.off()
  return(data)
}

lung1 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A1_KPx2_12min.rds', rctd = 'spatial_a1_RCTD_rna_reference.RData')
lung2 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A2_KPx2_6min.rds', rctd = 'spatial_a2_RCTD_rna_reference.RData')
lung3 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A3_KPx3_12min.rds', rctd = 'spatial_a3_RCTD_rna_reference.RData')
lung4 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A4_KPx3_6min.rds', rctd = 'spatial_a4_RCTD_rna_reference.RData')

cell_type_cor <- function(data1, data2){
  data <- merge(data1, data2, add.cell.ids = c('data1', 'data2'))
  matx0 <- t(data@assays$RCTD@data)
  correlation <- cor(matx0)
  print(Heatmap(correlation, show_row_dend = F, show_column_dend = F, show_column_names = F, ))
}

cell_type_cor_single <- function(data){
  matx0 <- t(data@assays$RCTD@data)
  correlation <- cor(matx0)
  print(Heatmap(correlation, show_row_dend = F, show_column_dend = F, show_column_names = F))
}

pdf('rna_reference_cell_type_cor_KPx2.pdf', width = 6, height = 4)
cell_type_cor(lung1, lung2)
dev.off()

pdf('rna_reference_cell_type_cor_KPx3.pdf', width = 6, height = 4)
cell_type_cor(lung3, lung4)
dev.off()

pdf('club_slide.pdf', width = 4, height = 4)
SpatialFeaturePlot(lung3, features = 'Club cells', pt.size.factor = 1.6, crop = TRUE)
dev.off()

pdf('Scgb1a1_slide.pdf', width = 4, height = 4)
DefaultAssay(lung3) <- 'SCT'
SpatialFeaturePlot(lung3, features = 'Scgb1a1', pt.size.factor = 1.6, crop = TRUE)
dev.off()

pdf('Scgb3a2_slide.pdf', width = 4, height = 4)
DefaultAssay(lung3) <- 'SCT'
SpatialFeaturePlot(lung3, features = 'Scgb3a2', pt.size.factor = 1.6, crop = TRUE)
dev.off()
