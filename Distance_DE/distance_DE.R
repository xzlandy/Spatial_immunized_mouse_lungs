library(Seurat)
library(openxlsx)
library(limma)
library(parallel)
library(reshape2)
library(splines)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(SplinesUtils)
library(circlize)

setwd('~/Box/RWorkSpace/Spatial_combine/Distance_DE/')

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
  cell <- data.matrix(mtx0[mtx0$Mgp < quantile(mtx0$Mgp, Mgp_cutoff) & mtx0[,cell] >= quantile(mtx0[,cell], cell_cutoff) , c("imagerow", "imagecol", cell)])
  cell.dist <- computeDist(cell, airCells)
  return(cell.dist)
}

cell_dist_true <- function(data, cell, cell_cutoff){
  mtx0 = data@images$slice1@coordinates
  mtx0 <- cbind(mtx0, t(data@assays$RCTD@data))
  mtx0 <- cbind(mtx0, data@meta.data)
  mtx0$Mgp <- data@assays$SCT@data['Mgp',]
  airCells <-  data.matrix(mtx0[mtx0$airway == 1, c("imagerow", "imagecol")])
  cell <- data.matrix(mtx0[mtx0$vessel == 0 & mtx0[,cell] >= quantile(mtx0[,cell], cell_cutoff) , c("imagerow", "imagecol", cell)])
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
      xlab('Distance to airway (micrometer)')
  }
  else if(type == 'immune'){
    cell <- cell[,-c(1:6, 17)]
    tmp <- melt(cell, id.vars = 'cell_dist')
    colnames(tmp)[2:3] <- c('Type', 'Value')
    
    ggplot(tmp, aes(x = cell_dist, y = Value, color = Type, linetype = Type))+
      geom_smooth(formula = y ~ ns(x, df = 3), se = T, alpha = 0.2)+
      theme_bw()+
      ylab('Proportions of cell types')+
      xlab('Distance to airway (micrometer)')
  }
}

plot_cell_dist_true_cutoff <- function(data, type, cutoff){
  mtx0 = data@images$slice1@coordinates
  mtx0 <- cbind(mtx0, t(data@assays$RCTD@data))
  mtx0 <- cbind(mtx0, data@meta.data)
  mtx0$Mgp <- data@assays$SCT@data['Mgp',]
  airCells <-  data.matrix(mtx0[mtx0$airway == 1, c("imagerow", "imagecol")])
  cell <- data.matrix(mtx0[mtx0$vessel == 0, c("imagerow", "imagecol", rownames(data@assays$RCTD@data))])
  cell.dist <- computeDist(cell, airCells)
  cell.dist <- cell.dist[cell.dist < cutoff]
  cell <- as.data.frame(cell)
  cell <- cell[names(cell.dist),]
  cell$cell_dist <- cell.dist
  
  if(type == 'all'){
    cell <- cell[,-c(1:2)]
    tmp <- melt(cell, id.vars = 'cell_dist')
    colnames(tmp)[2:3] <- c('Type', 'Value')
    
    ggplot(tmp, aes(x = cell_dist, y = Value, color = Type))+
      geom_smooth(formula = y ~ ns(x, df = 3), se = T, alpha = 0.2)+
      theme_bw()+
      ylab('Proportions of cell types')+
      xlab('Distance to airway (micrometer)')
  }
  else if(type == 'immune'){
    cell <- cell[,-c(1:6, 17)]
    tmp <- melt(cell, id.vars = 'cell_dist')
    colnames(tmp)[2:3] <- c('Type', 'Value')
    
    ggplot(tmp, aes(x = cell_dist, y = Value, color = Type, linetype = Type))+
      geom_smooth(formula = y ~ ns(x, df = 3), se = T, alpha = 0.2)+
      theme_bw()+
      ylab('Proportions of cell types')+
      xlab('Distance to airway (micrometer)')
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

distance_DE_field_true_spline <- function(raw_data){
  tmp <- cell_dist_true(raw_data, 'Th17', 0)
  
  raw_data$cell_dist <- NA
  raw_data@meta.data[names(tmp),]$cell_dist <- tmp
  
  coldata <- raw_data@meta.data
  countdata <- raw_data@assays$SCT@data
  coldata <- coldata[!is.na(coldata$cell_dist),]
  countdata <- countdata[,rownames(coldata)]
  proportion <- as.data.frame(t(raw_data@assays$RCTD@data))
  coldata <- cbind(coldata, t(raw_data@assays$RCTD@data)[rownames(coldata),])
  
  matrix <- model.matrix(~ns(cell_dist, df = 3), coldata)
  f1 <- lmFit(countdata, matrix)
  ef1 <- eBayes(f1)
  table1 <- topTable(ef1,2,number = Inf,adjust.method = "BH")
  table2 <- topTable(ef1,3,number = Inf,adjust.method = "BH")
  table2 <- table2[rownames(table1),]
  table3 <- topTable(ef1,4,number = Inf,adjust.method = "BH")
  table3 <- table3[rownames(table1),]
  table <- cbind(table1, table2)
  table <- cbind(table, table3)
  colnames(table) <- paste(colnames(table), rep(1:3, each = 6), sep = '_')
  table$gene <- rownames(table)
  table$P.Value_min <- pmin(table$P.Value_1, table$P.Value_2, table$P.Value_3)
  table$adj.P.Val_min <- pmin(table$adj.P.Val_1, table$adj.P.Val_2, table$adj.P.Val_3)
  # table$sig <- table$adj.P.Val_1 < 0.05 & table$adj.P.Val_2 < 0.05 & table$adj.P.Val_3 < 0.05
  table$sig <- table$adj.P.Val_min < 0.05
  table <- table[order(table$adj.P.Val_min),]
  print(table(table$sig))
  return(table)
}

distance_DE_field_true_spline_cutoff <- function(raw_data, cutoff){
  tmp <- cell_dist_true(raw_data, 'Th17', 0)
  
  raw_data$cell_dist <- NA
  raw_data@meta.data[names(tmp),]$cell_dist <- tmp
  raw_data <- raw_data[,raw_data$cell_dist < cutoff]
  
  coldata <- raw_data@meta.data
  countdata <- raw_data@assays$SCT@data
  coldata <- coldata[!is.na(coldata$cell_dist),]
  countdata <- countdata[,rownames(coldata)]
  proportion <- as.data.frame(t(raw_data@assays$RCTD@data))
  coldata <- cbind(coldata, t(raw_data@assays$RCTD@data)[rownames(coldata),])
  
  matrix <- model.matrix(~ns(cell_dist, df = 3), coldata)
  f1 <- lmFit(countdata, matrix)
  ef1 <- eBayes(f1)
  table1 <- topTable(ef1,2,number = Inf,adjust.method = "BH")
  table2 <- topTable(ef1,3,number = Inf,adjust.method = "BH")
  table2 <- table2[rownames(table1),]
  table3 <- topTable(ef1,4,number = Inf,adjust.method = "BH")
  table3 <- table3[rownames(table1),]
  table <- cbind(table1, table2)
  table <- cbind(table, table3)
  colnames(table) <- paste(colnames(table), rep(1:3, each = 6), sep = '_')
  table$gene <- rownames(table)
  table$P.Value_min <- pmin(table$P.Value_1, table$P.Value_2, table$P.Value_3)
  table$adj.P.Val_min <- pmin(table$adj.P.Val_1, table$adj.P.Val_2, table$adj.P.Val_3)
  # table$sig <- table$adj.P.Val_1 < 0.05 & table$adj.P.Val_2 < 0.05 & table$adj.P.Val_3 < 0.05
  table$sig <- table$adj.P.Val_min < 0.05
  table <- table[order(table$adj.P.Val_min),]
  print(table(table$sig))
  return(table)
}

# distance_DE_field_true_spline_nb <- function(raw_data, type){
#   tmp <- cell_dist_true(raw_data, 'Th17', 0)
#   
#   raw_data$cell_dist <- NA
#   raw_data@meta.data[names(tmp),]$cell_dist <- tmp
#   
#   coldata <- raw_data@meta.data
#   countdata <- raw_data@assays$Spatial@counts
#   coldata <- coldata[!is.na(coldata$cell_dist),]
#   countdata <- countdata[,rownames(coldata)]
#   proportion <- as.data.frame(t(raw_data@assays$RCTD@data))
#   coldata <- cbind(coldata, t(raw_data@assays$RCTD@data)[rownames(coldata),])
#   
#   zero_rate = apply(countdata, 1, function(x) mean(x == 0))
#   countdata = countdata[zero_rate != 1,]
#   table_init <- as.data.frame(matrix(NA, nrow = nrow(countdata), ncol = 4))
#   rownames(table_init) <- rownames(countdata)
#   table1=table2=table3=table4=table_init
#   library(parallel)
#   library(MASS)
#   library(pscl)
#   if(type == 'nb'){
#     for (i in 1:200){
#       if (i %% 100 == 0){
#         print(i)
#       }
#       l1 <- glm.nb(countdata[i,] ~ ns(coldata$cell_dist, df = 3))
#       coef <- summary(l1)$coefficients
#       table1[i,] <- coef[1,]
#       table2[i,] <- coef[2,]
#       table3[i,] <- coef[3,]
#       table4[i,] <- coef[4,]
#     }
#     table <- mclapply(1:100, function(i){
#       l1 <- glm.nb(countdata[i,] ~ ns(coldata$cell_dist, df = 3))
#       coef <- as.data.frame(summary(l1)$coefficients)
#       return(coef)
#     }, mc.cores = 2)
#   }
#   else if(type == 'zinb'){
#     mclapply(1:nrow(countdata), function(i){
#       l1 <- zinb(countdata[i,] ~ ns(coldata$cell_dist, df = 3))
#     })
#   }
# }

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

DE_combine <- function(marker1, marker2, name1, name2){
  marker1.sig <- marker1[marker1$sig == T,]
  marker2.sig <- marker2[marker2$sig == T,]
  merge <- inner_join(marker1.sig, marker2.sig, by = "gene", suffix = c(".A1",".A2"))
  merge$adj.P.Val_min <- pmin(merge$adj.P.Val_min.A1, merge$adj.P.Val_min.A2)
  merge <- merge[order(merge$adj.P.Val_min),]
  colnames(merge) <- gsub('.A1', name1, colnames(merge))
  colnames(merge) <- gsub('.A2', name2, colnames(merge))
  rownames(merge) <- merge$gene
  print(nrow(merge))
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

pdf('distance_cell_type_immune_A1.pdf', width = 6, height = 4)
plot_cell_dist_true(lung1, 'immune')
dev.off()

pdf('distance_cell_type_immune_A2.pdf', width = 6, height = 4)
plot_cell_dist_true(lung2, 'immune')
dev.off()

pdf('distance_cell_type_immune_A3.pdf', width = 6, height = 4)
plot_cell_dist_true(lung3, 'immune')
dev.off()

pdf('distance_cell_type_immune_A4.pdf', width = 6, height = 4)
plot_cell_dist_true(lung4, 'immune')
dev.off()

pdf('distance_cell_type_all_A1.pdf', width = 6, height = 4)
plot_cell_dist_true(lung1, 'all')
dev.off()

pdf('distance_cell_type_all_A2.pdf', width = 6, height = 4)
plot_cell_dist_true(lung2, 'all')
dev.off()

pdf('distance_cell_type_all_A3.pdf', width = 6, height = 4)
plot_cell_dist_true(lung3, 'all')
dev.off()

pdf('distance_cell_type_all_A4.pdf', width = 6, height = 4)
plot_cell_dist_true(lung4, 'all')
dev.off()

marker.1 <- distance_DE_field_true_spline(lung1)
marker.2 <- distance_DE_field_true_spline(lung2)
marker.3 <- distance_DE_field_true_spline(lung3)
marker.4 <- distance_DE_field_true_spline(lung4)

marker.KP2 <- DE_combine(marker.1, marker.2, '.A1', '.A2')
marker.KP3 <- DE_combine(marker.3, marker.4, '.A3', '.A4')

check_gene <- function(gene){
  SpatialFeaturePlot(merge.lung, paste0('sct_', gene)) / SpatialFeaturePlot(merge.lung, 'airway', pt.size.factor = 1.6, alpha = c(0,1))
}

check_gene_only <- function(gene){
  SpatialFeaturePlot(merge.lung, gene)
}

check_chemokine <- function(table){
  chemokine <- table[c(grep('^Cxcl', table$gene), grep('^Ccl', table$gene), grep('^Cx3cl', table$gene), grep('^Xcl', table$gene)),]
  chemokine <- chemokine[order(chemokine$adj.P.Val_min),]
  return(chemokine)
}

write.xlsx(marker.KP2, file = 'distance_DE_KP2.xlsx', rowNames = T)
write.xlsx(marker.KP3, file = 'distance_DE_KP3.xlsx', rowNames = T)

# distance_DE_field_true_spline <- function(raw_data1, raw_data2, genes){
#   raw_data <- raw_data1
#   tmp <- cell_dist_true(raw_data, 'Th17', 0)
#   
#   raw_data$cell_dist <- NA
#   raw_data@meta.data[names(tmp),]$cell_dist <- tmp
#   
#   coldata <- raw_data@meta.data
#   countdata <- raw_data@assays$SCT@data
#   coldata <- coldata[!is.na(coldata$cell_dist),]
#   countdata <- countdata[genes,rownames(coldata)]
#   proportion <- as.data.frame(t(raw_data@assays$RCTD@data))
#   coldata <- cbind(coldata, t(raw_data@assays$RCTD@data)[rownames(coldata),])
#   
#   cell_dist <- coldata$cell_dist
#   gene.list <- list()
#   predicted.table <- matrix(NA, nrow = nrow(countdata), ncol = length(seq(0, 1610, 1)))
#   rownames(predicted.table) <- rownames(countdata)
#   for (i in 1:nrow(countdata)){
#     f1 <- lm(countdata[i,] ~ ns(cell_dist, df = 3))
#     gene.list[[rownames(countdata)[i]]] <- f1
#     predicted.table[i,] <- predict(f1, data.frame(cell_dist = seq(0, 1610, 1)))
#   }
#   
#   Heatmap(predicted.table, cluster_columns = F, show_row_names = F)
#   return(table)
# }
View(check_chemokine(marker.KP3))
View(check_chemokine(marker.KP2))

write.xlsx(check_chemokine(marker.KP2), file = 'chemokine_KP2.xlsx', rowNames = T)
write.xlsx(check_chemokine(marker.KP3), file = 'chemokine_KP3.xlsx', rowNames = T)

check_gene('Cxcl17')
check_gene('Ccl20')
check_gene('Cxcl13')
check_gene('Ccl8')
check_gene('Cxcl12')
check_gene('Cxcl9')
check_gene('Ccl5')
check_gene('Ccl9')
check_gene('Ccl7')

DefaultAssay(merge.lung) <- 'SCT'
pdf('Ccl20.pdf', width = 16, height = 4)
check_gene_only('Ccl20')
dev.off()

merge <- merge(lung1, lung3)
names(merge@images) <- c("slice_A1","slice_A3")
pdf('Ccl20_2.pdf', width = 8, height = 4)
SpatialFeaturePlot(merge, 'Ccl20')
dev.off()

load('../Integration/T_cells/transferred_rna_reference.RData')

FeaturePlot(rna_reference, c('Cxcl17', 'Ccl20', 'Cxcl13', 'Ccl8', 'Cxcl12', 'Cxcl9', 'Ccl5', 'Ccl9', 'Ccl7'))

pdf('distance_cell_type_immune_A1_cutoff.pdf', width = 6, height = 4)
plot_cell_dist_true_cutoff(lung1, 'immune', 1000)
dev.off()

pdf('distance_cell_type_immune_A2_cutoff.pdf', width = 6, height = 4)
plot_cell_dist_true_cutoff(lung2, 'immune', 1000)
dev.off()

pdf('distance_cell_type_immune_A3_cutoff.pdf', width = 6, height = 4)
plot_cell_dist_true_cutoff(lung3, 'immune', 1000)
dev.off()

pdf('distance_cell_type_immune_A4_cutoff.pdf', width = 6, height = 4)
plot_cell_dist_true_cutoff(lung4, 'immune', 1000)
dev.off()

pdf('distance_cell_type_all_A1_cutoff.pdf', width = 6, height = 4)
plot_cell_dist_true_cutoff(lung1, 'all', 1000)
dev.off()

pdf('distance_cell_type_all_A2_cutoff.pdf', width = 6, height = 4)
plot_cell_dist_true_cutoff(lung2, 'all', 1000)
dev.off()

pdf('distance_cell_type_all_A3_cutoff.pdf', width = 6, height = 4)
plot_cell_dist_true_cutoff(lung3, 'all', 1000)
dev.off()

pdf('distance_cell_type_all_A4_cutoff.pdf', width = 6, height = 4)
plot_cell_dist_true_cutoff(lung4, 'all', 1000)
dev.off()

marker.1_cutoff <- distance_DE_field_true_spline_cutoff(lung1, 1000)
marker.2_cutoff <- distance_DE_field_true_spline_cutoff(lung2, 1000)
marker.3_cutoff <- distance_DE_field_true_spline_cutoff(lung3, 1000)
marker.4_cutoff <- distance_DE_field_true_spline_cutoff(lung4, 1000)

marker.KP2_cutoff <- DE_combine(marker.1_cutoff, marker.2_cutoff, '.A1', '.A2')
marker.KP3_cutoff <- DE_combine(marker.3_cutoff, marker.4_cutoff, '.A3', '.A4')

write.xlsx(marker.KP2_cutoff, file = 'distance_DE_KP2_cutoff.xlsx', rowNames = T)
write.xlsx(marker.KP3_cutoff, file = 'distance_DE_KP3_cutoff.xlsx', rowNames = T)

write.xlsx(check_chemokine(marker.KP2_cutoff), file = 'chemokine_KP2_cutoff.xlsx', rowNames = T)
write.xlsx(check_chemokine(marker.KP3_cutoff), file = 'chemokine_KP3_cutoff.xlsx', rowNames = T)
