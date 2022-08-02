library(Seurat)
library(openxlsx)
library(limma)
library(parallel)
library(reshape2)
library(splines)
library(ggplot2)

setwd('~/Box/RWorkSpace/Spatial_combine/DE_analysis/')

load('../Spatial/Rfiles/lung_spatial.RData')
load('../Spatial/Rfiles/merge_lung_spatial.RData')

computeDist = function(subCell, airCells) {
  nAir = nrow(airCells)
  helper = matrix(rep(1,nrow(subCell)), ncol=1)
  
  outMtx = matrix(nrow = nrow(subCell), ncol = nAir)
  rownames(outMtx) = rownames(subCell)
  colnames(outMtx) = rownames(airCells)
  
  for (i in 1:nAir) {
    temp0 = subCell[,c(1,2)] - helper %*% airCells[i,]
    outMtx[,i] = sqrt(temp0[,1]^2 + temp0[,2]^2)
  }
  min.dis = apply(outMtx,1,min)
  return(min.dis)
}

cell_dist <- function(data, club_cutoff, Mgp_cutoff, cell, cell_cutoff){
  mtx0 = data@images$slice1@coordinates
  mtx0 <- cbind(mtx0, t(data@assays$RCTD@data))
  mtx0 <- cbind(mtx0, data@meta.data)
  mtx0$Mgp <- data@assays$SCT@data['Mgp',]
  airCells <-  data.matrix(mtx0[mtx0$`Club cells` > quantile(mtx0$`Club cells`, club_cutoff), c("imagerow", "imagecol")])
  cell <- data.matrix(mtx0[mtx0$Mgp < quantile(mtx0$Mgp, Mgp_cutoff) & mtx0[,cell] > quantile(mtx0[,cell], cell_cutoff) , c("imagerow", "imagecol", cell)])
  cell.dist <- computeDist(cell, airCells)
  return(cell.dist)
}

cell_dist_true <- function(data, cell, cell_cutoff){
  mtx0 = data@images$slice1@coordinates
  mtx0 <- cbind(mtx0, t(data@assays$RCTD@data))
  mtx0 <- cbind(mtx0, data@meta.data)
  mtx0$Mgp <- data@assays$SCT@data['Mgp',]
  airCells <-  data.matrix(mtx0[mtx0$airway == 1, c("imagerow", "imagecol")])
  cell <- data.matrix(mtx0[mtx0$vessel == 0 & mtx0[,cell] > quantile(mtx0[,cell], cell_cutoff) , c("imagerow", "imagecol", cell)])
  cell.dist <- computeDist(cell, airCells)
  return(cell.dist)
}

plot_cell_dist_true <- function(data, type){
  mtx0 = data@images$slice1@coordinates
  mtx0 <- cbind(mtx0, t(data@assays$RCTD@data))
  mtx0 <- cbind(mtx0, data@meta.data)
  mtx0$Mgp <- data@assays$SCT@data['Mgp',]
  airCells <-  data.matrix(mtx0[mtx0$airway == 1, c("imagerow", "imagecol")])
  cell <- data.matrix(mtx0[mtx0$vessel == 0, c("imagerow", "imagecol", rownames(data@assays$RCTD@data))])
  cell.dist <- computeDist(cell, airCells)
  cell <- as.data.frame(cell)
  cell$cell_dist <- cell.dist
  
  if(type == 'all'){
    cell <- cell[,-c(1:2)]
    tmp <- melt(cell, id.vars = 'cell_dist')
    colnames(tmp)[2:3] <- c('Type', 'Value')
    
    ggplot(tmp, aes(x = cell_dist, y = Value, color = Type))+
      geom_smooth(formula = y ~ ns(x, df = 3), se = T, alpha = 0.2)+
      theme_bw()+
      ylab('Proportions of cell types')+
      xlab('Distance to airway')
  }
  else if(type == 'immune'){
    cell <- cell[,-c(1:6, 17)]
    tmp <- melt(cell, id.vars = 'cell_dist')
    colnames(tmp)[2:3] <- c('Type', 'Value')
    
    ggplot(tmp, aes(x = cell_dist, y = Value, color = Type, linetype = Type))+
      geom_smooth(formula = y ~ ns(x, df = 3), se = T, alpha = 0.2)+
      theme_bw()+
      ylab('Proportions of cell types')+
      xlab('Distance to airway')
  }
}

distance_DE_field <- function(raw_data){
  tmp <- cell_dist(raw_data, 0.90, 0.90, 'Th17', 0)
  
  raw_data$cell_dist <- NA
  raw_data@meta.data[names(tmp),]$cell_dist <- tmp
  
  coldata <- raw_data@meta.data
  countdata <- raw_data@assays$SCT@data
  coldata <- coldata[!is.na(coldata$cell_dist),]
  countdata <- countdata[,rownames(coldata)]
  proportion <- as.data.frame(t(raw_data@assays$RCTD@data))
  coldata <- cbind(coldata, t(raw_data@assays$RCTD@data)[rownames(coldata),])
  
  matrix <- model.matrix(~cell_dist, coldata)
  f1 <- lmFit(countdata, matrix)
  ef1 <- eBayes(f1)
  table <- topTable(ef1,2,number = Inf,adjust.method = "fdr")
  table$gene <- rownames(table)
  return(table)
}

distance_DE_field_true <- function(raw_data){
  tmp <- cell_dist_true(raw_data, 'Th17', 0)
  
  raw_data$cell_dist <- NA
  raw_data@meta.data[names(tmp),]$cell_dist <- tmp
  
  coldata <- raw_data@meta.data
  countdata <- raw_data@assays$SCT@data
  coldata <- coldata[!is.na(coldata$cell_dist),]
  countdata <- countdata[,rownames(coldata)]
  proportion <- as.data.frame(t(raw_data@assays$RCTD@data))
  coldata <- cbind(coldata, t(raw_data@assays$RCTD@data)[rownames(coldata),])
  
  matrix <- model.matrix(~cell_dist, coldata)
  f1 <- lmFit(countdata, matrix)
  ef1 <- eBayes(f1)
  table <- topTable(ef1,2,number = Inf,adjust.method = "fdr")
  table$gene <- rownames(table)
  return(table)
}

distance_DE <- function(raw_data, cell, cell_cutoff){
  tmp <- cell_dist(raw_data, 0.90, 0.90, cell, cell_cutoff)
  tmp <- tmp[tmp > 0]
  
  raw_data$cell_dist <- NA
  raw_data@meta.data[names(tmp),]$cell_dist <- tmp
  
  coldata <- raw_data@meta.data
  countdata <- raw_data@assays$SCT@data
  coldata <- coldata[!is.na(coldata$cell_dist),]
  countdata <- countdata[,rownames(coldata)]
  proportion <- as.data.frame(t(raw_data@assays$RCTD@data))
  coldata <- cbind(coldata, t(raw_data@assays$RCTD@data)[rownames(coldata),])
  
  matrix <- model.matrix(~cell_dist, coldata)
  f1 <- lmFit(countdata, matrix)
  ef1 <- eBayes(f1)
  table <- topTable(ef1,2,number = Inf,adjust.method = "fdr")
  table$gene <- rownames(table)
  return(table)
}

distance_DE_true <- function(raw_data, cell, cell_cutoff){
  tmp <- cell_dist_true(raw_data, cell, cell_cutoff)
  tmp <- tmp[tmp > 0]
  
  raw_data$cell_dist <- NA
  raw_data@meta.data[names(tmp),]$cell_dist <- tmp
  
  coldata <- raw_data@meta.data
  countdata <- raw_data@assays$SCT@data
  coldata <- coldata[!is.na(coldata$cell_dist),]
  countdata <- countdata[,rownames(coldata)]
  proportion <- as.data.frame(t(raw_data@assays$RCTD@data))
  coldata <- cbind(coldata, t(raw_data@assays$RCTD@data)[rownames(coldata),])
  
  matrix <- model.matrix(~cell_dist, coldata)
  f1 <- lmFit(countdata, matrix)
  ef1 <- eBayes(f1)
  table <- topTable(ef1,2,number = Inf,adjust.method = "fdr")
  table$gene <- rownames(table)
  return(table)
}

DE_meta <- function(marker1, marker2, name1, name2){
  merge <- merge(marker1, marker2, by = "gene", suffixes = c(".A1",".A2"))
  merge <- merge[,-c(7,13)]
  w.A1 <- 1/(merge$logFC.A1/merge$t.A1)^2
  w.A2 <- 1/(merge$logFC.A2/merge$t.A2)^2
  SE.meta <- sqrt(1/(w.A1 + w.A2))
  merge$logFC.meta <- (merge$logFC.A1 * w.A1 + merge$logFC.A2 * w.A2)/(w.A1 + w.A2)
  z.meta <- merge$logFC.meta/SE.meta
  merge$P.Value.meta <- 2 * pnorm(-abs(z.meta))
  merge$adj.P.Val.meta <- p.adjust(merge$P.Value.meta, method="fdr")
  merge <- merge[order(merge$P.Value.meta),]
  merge <- merge[,-c(4,9)]
  merge$both_sig <- merge$adj.P.Val.A1 < 0.05 & merge$adj.P.Val.A2 < 0.05
  merge$same_dir <- merge$logFC.A1 * merge$logFC.A2 > 0
  colnames(merge) <- gsub('.A1', name1, colnames(merge))
  colnames(merge) <- gsub('.A2', name2, colnames(merge))
  rownames(merge) <- 1:nrow(merge)
  print(nrow(merge[merge$adj.P.Val.meta < 0.05,]))
  print(table(merge$both_sig))
  return(merge)
}

check_contour <- function(type){
  p1 <- SpatialFeaturePlot(lung1, type, pt.size.factor = 1.6, alpha = c(0,1))
  p2 <- SpatialFeaturePlot(lung2, type, pt.size.factor = 1.6, alpha = c(0,1))
  p3 <- SpatialFeaturePlot(lung3, type, pt.size.factor = 1.6, alpha = c(0,1))
  p4 <- SpatialFeaturePlot(lung4, type, pt.size.factor = 1.6, alpha = c(0,1))
  plot <- (p1 | p2 | p3 | p4)
  return(plot)
}

DimPlot(merge.lung, reduction = "umap", group.by = c("ident", "orig.ident"))
FeaturePlot(merge.lung, 'rctd_Club cells') | DimPlot(merge.lung, reduction = "umap", group.by = c("orig.ident"))

check_gene <- function(gene){
  SpatialFeaturePlot(merge.lung, paste0('sct_', gene)) / SpatialFeaturePlot(merge.lung, 'airway', pt.size.factor = 1.6, alpha = c(0,1))
}

check_gene_only <- function(gene){
  SpatialFeaturePlot(merge.lung, gene)
}

load('../Integration/T_cells/transferred_rna_reference.RData')

Idents(merge.lung) <- 'orig.ident'
merge.lung <- RenameIdents(merge.lung,
                           'A1' = 'KPx2',
                           'A2' = 'KPx2',
                           'A3' = 'KPx3',
                           'A4' = 'KPx3')
merge.lung$conditions <- Idents(merge.lung)

merge.lung$airway_type <- 'Other'
merge.lung@meta.data[merge.lung$airway == 1 & merge.lung$conditions == 'KPx3',]$airway_type <- 'KPx3_airway'
merge.lung@meta.data[merge.lung$airway == 1 & merge.lung$conditions == 'KPx2',]$airway_type <- 'KPx2_airway'

Idents(merge.lung) <- 'airway_type'
DefaultAssay(merge.lung) <- 'SCT'
marker.airway_wilcox <- FindMarkers(merge.lung, ident.1 = 'KPx3_airway', ident.2 = 'KPx2_airway', logfc.threshold = 0, min.pct = 0)
marker.airway_wilcox$gene <- rownames(marker.airway_wilcox)
marker.airway_mast <- FindMarkers(merge.lung, ident.1 = 'KPx3_airway', ident.2 = 'KPx2_airway', test.use = "MAST", logfc.threshold = 0, min.pct = 0)
marker.airway_mast$gene <- rownames(marker.airway_mast)

pdf('AW112010.pdf', width = 16, height = 4)
check_gene_only('AW112010')
dev.off()

pdf('Retnla.pdf', width = 16, height = 4)
check_gene_only('Retnla')
dev.off()

pdf('Cbr2.pdf', width = 16, height = 4)
check_gene_only('Cbr2')
dev.off()

pdf('Cyp2f2.pdf', width = 16, height = 4)
check_gene_only('Cyp2f2')
dev.off()

write.xlsx(marker.airway_wilcox, file = 'airway_DE_wilcox.xlsx', rowNames = T)
write.xlsx(marker.airway_mast, file = 'airway_DE_mast.xlsx', rowNames = T)

define_cell <- function(raw_data){
  proportion <- as.data.frame(t(raw_data@assays$RCTD@data))
  colnames(proportion) <- gsub(' ', '_', colnames(proportion))
  cells <- round(colMeans(proportion)*nrow(proportion))
  
  tmp <- proportion %>% top_n(cells[1], get(names(cells[1])))
  raw_data$Alveolar_epithelial_cells <- 0
  raw_data$Alveolar_epithelial_cells[rownames(tmp)] <- 1
  
  tmp <- proportion %>% top_n(cells[2], get(names(cells[2])))
  raw_data$Club_cells <- 0
  raw_data$Club_cells[rownames(tmp)] <- 1
  
  tmp <- proportion %>% top_n(cells[5], get(names(cells[5])))
  raw_data$Monocytes <- 0
  raw_data$Monocytes[rownames(tmp)] <- 1
  
  tmp <- proportion %>% top_n(cells[6], get(names(cells[6])))
  raw_data$Macrophages <- 0
  raw_data$Macrophages[rownames(tmp)] <- 1
  
  tmp <- proportion %>% top_n(cells[7], get(names(cells[7])))
  raw_data$Dendritic_cells <- 0
  raw_data$Dendritic_cells[rownames(tmp)] <- 1
  
  tmp <- proportion %>% top_n(cells[8], get(names(cells[8])))
  raw_data$Neutrophils <- 0
  raw_data$Neutrophils[rownames(tmp)] <- 1
  
  tmp <- proportion %>% top_n(cells[9], get(names(cells[9])))
  raw_data$B_cells <- 0
  raw_data$B_cells[rownames(tmp)] <- 1
  
  tmp <- proportion %>% top_n(cells[10], get(names(cells[10])))
  raw_data$Th17 <- 0
  raw_data$Th17[rownames(tmp)] <- 1
  
  tmp <- proportion %>% top_n(cells[11], get(names(cells[11])))
  raw_data$Th1 <- 0
  raw_data$Th1[rownames(tmp)] <- 1
  
  tmp <- proportion %>% top_n(cells[12], get(names(cells[12])))
  raw_data$Other_T_cells_1 <- 0
  raw_data$Other_T_cells_1[rownames(tmp)] <- 1
  
  tmp <- proportion %>% top_n(cells[13], get(names(cells[13])))
  raw_data$Other_T_cells_2 <- 0
  raw_data$Other_T_cells_2[rownames(tmp)] <- 1
  
  tmp <- proportion %>% top_n(cells[14], get(names(cells[14])))
  raw_data$NK_cells <- 0
  raw_data$NK_cells[rownames(tmp)] <- 1
  
  raw_data$cells_no <- (raw_data$Alveolar_epithelial_cells + raw_data$Club_cells + raw_data$Monocytes + raw_data$Macrophages 
                        + raw_data$Dendritic_cells + raw_data$Neutrophils + raw_data$B_cells 
                        + raw_data$Th17 + raw_data$Th1 + raw_data$Other_T_cells_1
                        + raw_data$Other_T_cells_2 + raw_data$NK_cells)
  
  raw_data$cells <- raw_data$cells_no
  raw_data$cells[raw_data$cells > 1] <- 0
  
  raw_data$cell_types <- NA
  raw_data$cell_types[raw_data$Alveolar_epithelial_cells == 1 & raw_data$cells == 1] <- 'Alveolar epithelial cells'
  raw_data$cell_types[raw_data$Club_cells == 1 & raw_data$cells == 1] <- 'Club cells'
  raw_data$cell_types[raw_data$Monocytes == 1 & raw_data$cells == 1] <- 'Monocytes'
  raw_data$cell_types[raw_data$Macrophages == 1 & raw_data$cells == 1] <- 'Macrophages'
  raw_data$cell_types[raw_data$Dendritic_cells == 1 & raw_data$cells == 1] <- 'Dendritic cells'
  raw_data$cell_types[raw_data$Neutrophils == 1 & raw_data$cells == 1] <- 'Neutrophils'
  raw_data$cell_types[raw_data$B_cells == 1 & raw_data$cells == 1] <- 'B cells'
  raw_data$cell_types[raw_data$Th17 == 1 & raw_data$cells == 1] <- 'Th17'
  raw_data$cell_types[raw_data$Th1 == 1 & raw_data$cells == 1] <- 'Th1'
  raw_data$cell_types[raw_data$Other_T_cells_1 == 1 & raw_data$cells == 1] <- 'Other T cells 1'
  raw_data$cell_types[raw_data$Other_T_cells_2 == 1 & raw_data$cells == 1] <- 'Other T cells 2'
  raw_data$cell_types[raw_data$NK_cells == 1 & raw_data$cells == 1] <- 'NK cells'
  
  Idents(raw_data) <- 'cell_types'
  return(raw_data)
}

lung1 <- define_cell(lung1)
lung2 <- define_cell(lung2)
lung3 <- define_cell(lung3)
lung4 <- define_cell(lung4)

merge.lung$cell_types <- NA
merge.lung@meta.data[paste('A1', colnames(lung1), sep = '_'),]$cell_types <- lung1$cell_types
merge.lung@meta.data[paste('A2', colnames(lung2), sep = '_'),]$cell_types <- lung2$cell_types
merge.lung@meta.data[paste('A3', colnames(lung3), sep = '_'),]$cell_types <- lung3$cell_types
merge.lung@meta.data[paste('A4', colnames(lung4), sep = '_'),]$cell_types <- lung4$cell_types
merge.lung$cell_types <- factor(merge.lung$cell_types, levels = c('Alveolar epithelial cells', 'Club cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'Th17', 'Th1', 'Other T cells 1', 'Other T cells 2', 'NK cells'))

merge.lung$mouse <- NA
merge.lung$mouse[merge.lung$conditions == 'KPx2'] <- 'Immunized Mouse'
merge.lung$mouse[merge.lung$conditions == 'KPx3'] <- 'Re-challenged Mouse'

Idents(merge.lung) <- 'orig.ident'
pdf('UMAP_batch.pdf', width = 8, height = 4)
DimPlot(merge.lung, split.by = 'mouse')
dev.off()

pdf('UMAP_airway.pdf', width = 8, height = 4)
FeaturePlot(merge.lung, 'airway', split.by = 'mouse')
dev.off()

pdf('UMAP_airway_slice.pdf', width = 16, height = 4)
FeaturePlot(merge.lung, 'airway', split.by = 'orig.ident')
dev.off()

Idents(merge.lung) <- 'cell_types'
pdf('cell_types.pdf', width = 32, height = 8)
SpatialDimPlot(merge.lung[,!is.na(merge.lung$cell_types)])
dev.off()

pdf('cell_types_new.pdf', width = 16, height = 4)
SpatialDimPlot(merge.lung[,!is.na(merge.lung$cell_types)])
dev.off()

pdf('UMAP_cell_types.pdf', width = 8, height = 4)
DimPlot(merge.lung[,!is.na(merge.lung$cell_types)], split.by = 'mouse')
dev.off()

pdf('UMAP_cell_types_slice.pdf', width = 16, height = 4)
DimPlot(merge.lung[,!is.na(merge.lung$cell_types)], split.by = 'orig.ident')
dev.off()

Idents(merge.lung) <- 'orig.ident'
pdf('UMAP_slice_cell_types.pdf', width = 48, height = 4)
DimPlot(merge.lung[,!is.na(merge.lung$cell_types)], split.by = 'cell_types')
dev.off()

Idents(merge.lung) <- 'cell_types'
pdf('dotplot_markers_final.pdf', width = 6, height = 12)
p <- DotPlot(merge.lung[,!is.na(merge.lung$cell_types)], features = c('Sftpa1', 'Sftpb', 'Sftpc', 'Sftpd',
                                         'Scgb1a1', 'Muc5b', 'Scgb3a1', 'Scgb3a2',
                                         'Cd14', 'Itgam',
                                         'Itgax', 'Cd68', 'Mrc1', 'Marco',
                                         'Aif1', 'H2-DMb1', 'H2-Eb1', 'H2-Aa',
                                         'Gsr', 'Pglyrp1', 'S100a8', 'Ly6g',
                                         'Ms4a1', 'Cd79a', 'Igkc', 'Cd19',
                                         'Cd3d', 'Cd3e', 'Cd3g',
                                         'Rora', 'Rorc', 'Ccr6', 'Ccr4',
                                         'Ifng', 'Tbx21', 'Cxcr5', 'Ccr5', 'Il12rb2',
                                         'Ncr1', 'Nkg7'))
p + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 0), axis.text.y = element_text(face = 'italic', hjust = 0), legend.position = 'bottom', legend.direction = 'vertical') + xlab('') + ylab('')
dev.off()

DE_cell_types <- function(cell_types, name){
  merge.lung$test_types <- 'Other'
  merge.lung@meta.data[merge.lung$cell_types %in% cell_types & merge.lung$conditions == 'KPx3',]$test_types <- 'KPx3_test_types'
  merge.lung@meta.data[merge.lung$cell_types %in% cell_types & merge.lung$conditions == 'KPx2',]$test_types <- 'KPx2_test_types'
  
  Idents(merge.lung) <- 'test_types'
  DefaultAssay(merge.lung) <- 'SCT'
  marker.wilcox <- FindMarkers(merge.lung, ident.1 = 'KPx3_test_types', ident.2 = 'KPx2_test_types', logfc.threshold = 0, min.pct = 0)
  marker.wilcox$gene <- rownames(marker.wilcox)
  marker.mast <- FindMarkers(merge.lung, ident.1 = 'KPx3_test_types', ident.2 = 'KPx2_test_types', test.use = "MAST", logfc.threshold = 0, min.pct = 0)
  marker.mast$gene <- rownames(marker.mast)
  
  write.xlsx(marker.wilcox, file = paste(name, 'DE_wilcox.xlsx', sep = '_'), rowNames = T)
  write.xlsx(marker.mast, file = paste(name, 'DE_mast.xlsx', sep = '_'), rowNames = T)
}

Idents(merge.lung) <- 'conditions'
DefaultAssay(merge.lung) <- 'SCT'
marker.wilcox <- FindMarkers(merge.lung, ident.1 = 'KPx3', ident.2 = 'KPx2', logfc.threshold = 0, min.pct = 0)
marker.wilcox$gene <- rownames(marker.wilcox)
marker.mast <- FindMarkers(merge.lung, ident.1 = 'KPx3', ident.2 = 'KPx2', test.use = "MAST", logfc.threshold = 0, min.pct = 0)
marker.mast$gene <- rownames(marker.mast)

write.xlsx(marker.wilcox, file = 'conditions_DE_wilcox.xlsx', rowNames = T)
write.xlsx(marker.mast, file = 'conditions_DE_mast.xlsx', rowNames = T)

DE_cell_types(c('Th1', 'Th17', 'Other T cells 1', 'Other T cells 2'), 'T_cells')
DE_cell_types('B cells', 'B_cells')
DE_cell_types('Dendritic cells', 'Dendritic_cells')
DE_cell_types('Macrophages', 'Macrophages')
DE_cell_types('Monocytes', 'Monocytes')
DE_cell_types('Neutrophils', 'Neutrophils')
DE_cell_types('NK cells', 'NK_cells')
DE_cell_types('Club cells', 'Club_cells')
DE_cell_types('Th1', 'Th1')
DE_cell_types('Th17', 'Th17')
DE_cell_types('Other T cells 1', 'Other_T_cells_1')
DE_cell_types('Other T cells 2', 'Other_T_cells_2')
