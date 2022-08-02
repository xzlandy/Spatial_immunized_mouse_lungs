library(RCTD)
library(Matrix)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(ggpubr)
library(ggplot2)

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
  pdf(gsub('.RData', '.pdf', rctd), width = 20, height = 12)
  print(SpatialFeaturePlot(data, features = rownames(data), pt.size.factor = 1.6, ncol = 5, crop = TRUE))
  dev.off()
  return(data)
}

lung1 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A1_KPx2_12min.rds', rctd = 'spatial_a1_RCTD_transferred_rna_reference.RData')
lung2 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A2_KPx2_6min.rds', rctd = 'spatial_a2_RCTD_transferred_rna_reference.RData')
lung3 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A3_KPx3_12min.rds', rctd = 'spatial_a3_RCTD_transferred_rna_reference.RData')
lung4 <- read_spatial_RCTD(spatial = '../Spatial/Rfiles/A4_KPx3_6min.rds', rctd = 'spatial_a4_RCTD_transferred_rna_reference.RData')

cell_type_cor <- function(data1, data2){
  data <- merge(data1, data2, add.cell.ids = c('data1', 'data2'))
  matx0 <- t(data@assays$RCTD@data)
  correlation <- cor(matx0)
  print(Heatmap(correlation, show_row_dend = F, show_column_dend = F, show_column_names = F, heatmap_legend_param = list(title = "Correlation")))
}

cell_type_cor_single <- function(data){
  matx0 <- t(data@assays$RCTD@data)
  correlation <- cor(matx0)
  print(Heatmap(correlation, show_row_dend = F, show_column_dend = F, show_column_names = F, heatmap_legend_param = list(title = "Correlation")))
}

pdf('transferred_rna_reference_cell_type_cor_KPx2.pdf', width = 6, height = 4)
cell_type_cor(lung1, lung2)
dev.off()

pdf('transferred_rna_reference_cell_type_cor_KPx3.pdf', width = 6, height = 4)
cell_type_cor(lung3, lung4)
dev.off()

pdf('transferred_rna_reference_cell_type_cor_A1.pdf', width = 6, height = 4)
cell_type_cor_single(lung1)
dev.off()

pdf('transferred_rna_reference_cell_type_cor_A2.pdf', width = 6, height = 4)
cell_type_cor_single(lung2)
dev.off()

pdf('transferred_rna_reference_cell_type_cor_A3.pdf', width = 6, height = 4)
cell_type_cor_single(lung3)
dev.off()

pdf('transferred_rna_reference_cell_type_cor_A4.pdf', width = 6, height = 4)
cell_type_cor_single(lung4)
dev.off()

library(reshape2)
proportion <- function(data, name){
  matx0 <- as.data.frame(t(data@assays$RCTD@data))
  matx0$Slice <- name
  table <- melt(matx0, id.vars = 'Slice')
  colnames(table)[2:3] <- c('Cell', 'Proportion')
  return(table)
}

data.plot <- proportion(lung1, 'A1')
data.plot <- rbind(data.plot, proportion(lung2, 'A2'))
data.plot <- rbind(data.plot, proportion(lung3, 'A3'))
data.plot <- rbind(data.plot, proportion(lung4, 'A4'))
data.plot$Proportion <- data.plot$Proportion*100

data.plot$Condition <- NA
data.plot$Condition[data.plot$Slice %in% c('A1', 'A2')] <- 'KPx2'
data.plot$Condition[data.plot$Slice %in% c('A3', 'A4')] <- 'KPx3'

pdf('transferred_rna_reference_proportion_pvalue.pdf', width = 6, height = 12)
p <- ggboxplot(data.plot, x = "Condition", y = "Proportion",
               fill = "Slice", palette = 'npg', outlier.shape = NA)+
  ylim(0,60)+
  coord_flip()+
  rremove('y.text')
facet(p, facet.by = 'Cell', ncol = 1, strip.position = 'left')+
  stat_compare_means(method = 't.test', label = 'p.signif')+
  theme(strip.text.y.left = element_text(angle = 0, size = 14)) + ylab('Proportion (%)') + xlab('')
dev.off()

library(ggpubr)
pdf('transferred_rna_reference_proportion.pdf', width = 6, height = 6)
ggboxplot(data.plot, x = "Cell", y = "Proportion",
               fill = "Slice", palette = 'npg', outlier.shape = NA)+
  ylim(0,60)+
  coord_flip() + theme(axis.text.y = element_text(hjust = 0)) + ylab('Proportion (%)') + xlab('')
dev.off()
