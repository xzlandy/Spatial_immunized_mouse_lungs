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
  pdf(gsub('.RData', '.pdf', rctd), width = 20, height = 16)
  print(SpatialFeaturePlot(data, features = rownames(data), pt.size.factor = 1.6, ncol = 5, crop = TRUE))
  dev.off()
  return(data)
}

lung1 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A1_KPx2_12min.rds', rctd = 'spatial_a1_RCTD_cohen_lung.RData')
lung2 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A2_KPx2_6min.rds', rctd = 'spatial_a2_RCTD_cohen_lung.RData')
lung3 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A3_KPx3_12min.rds', rctd = 'spatial_a3_RCTD_cohen_lung.RData')
lung4 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A4_KPx3_6min.rds', rctd = 'spatial_a4_RCTD_cohen_lung.RData')

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

pdf('cohen_lung_cell_type_cor_KPx2.pdf', width = 6, height = 4)
cell_type_cor(lung1, lung2)
dev.off()

pdf('cohen_lung_cell_type_cor_KPx3.pdf', width = 6, height = 4)
cell_type_cor(lung3, lung4)
dev.off()

read_spatial_RCTD_to_compare <- function(spatial, rna_reference, cohen_lung){
  data <- readRDS(spatial)
  load(rna_reference)
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
  
  load(cohen_lung)
  results <- myRCTD@results
  norm_weights = as.matrix(sweep(results$weights, 1, rowSums(results$weights), '/'))
  norm_weights <- t(norm_weights)
  ncol <- length(colnames(data)[!(colnames(data) %in% colnames(norm_weights))])
  if (ncol > 0) {
    tmp <- matrix(0, nrow = nrow(norm_weights), ncol = ncol)
    colnames(tmp) <- colnames(data)[!(colnames(data) %in% colnames(norm_weights))]
    norm_weights <- cbind(norm_weights, tmp)[,colnames(data)]
  }
  data[['RCTD_Cohen']] <- CreateAssayObject(data = norm_weights)
  return(data)
}

################
lung1 <- read_spatial_RCTD_to_compare(spatial = '../Spatial/Rfiles/A1_KPx2_12min.rds', rna_reference = 'spatial_a1_RCTD_rna_reference.RData', cohen_lung = 'spatial_a1_RCTD_cohen_lung.RData')
lung2 <- read_spatial_RCTD_to_compare(spatial = '../Spatial/Rfiles/A2_KPx2_6min.rds', rna_reference = 'spatial_a2_RCTD_rna_reference.RData', cohen_lung = 'spatial_a2_RCTD_cohen_lung.RData')
lung3 <- read_spatial_RCTD_to_compare(spatial = '../Spatial/Rfiles/A3_KPx3_12min.rds', rna_reference = 'spatial_a3_RCTD_rna_reference.RData', cohen_lung = 'spatial_a3_RCTD_cohen_lung.RData')
lung4 <- read_spatial_RCTD_to_compare(spatial = '../Spatial/Rfiles/A4_KPx3_6min.rds', rna_reference = 'spatial_a4_RCTD_rna_reference.RData', cohen_lung = 'spatial_a4_RCTD_cohen_lung.RData')

deconvolution_to_compare <- function(data){
  matx1 <- t(data@assays$RCTD@data)
  matx2 <- t(data@assays$RCTD_Cohen@data)
  correlation <- cor(matx1, matx2)
  print(Heatmap(correlation, show_row_dend = F, show_column_dend = F, column_title = 'Public scRNA-seq Reference', row_title = 'In-house scRNA-seq Reference', heatmap_legend_param = list(title = "Correlation")))
}

pdf('deconvolution_to_compare_A1.pdf', width = 6, height = 4)
deconvolution_to_compare(lung1)
dev.off()

pdf('deconvolution_to_compare_A2.pdf', width = 6, height = 4)
deconvolution_to_compare(lung2)
dev.off()

pdf('deconvolution_to_compare_A3.pdf', width = 6, height = 4)
deconvolution_to_compare(lung3)
dev.off()

pdf('deconvolution_to_compare_A4.pdf', width = 6, height = 4)
deconvolution_to_compare(lung4)
dev.off()

pdf('tcell_rna_reference.pdf', width = 4, height = 4)
DefaultAssay(lung3) <- 'RCTD'
SpatialFeaturePlot(lung3, features = 'T cells', pt.size.factor = 1.6, crop = TRUE)
dev.off()

pdf('tcell_cohen.pdf', width = 4, height = 4)
DefaultAssay(lung3) <- 'RCTD_Cohen'
SpatialFeaturePlot(lung3, features = 'T cell', pt.size.factor = 1.6, crop = TRUE)
dev.off()