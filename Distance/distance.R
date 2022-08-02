library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(parallel)
library(ggpubr)

setwd('~/Box/RWorkSpace/Spatial_combine/Distance/')

# Function to compute distance
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

th17_th1_dist <- function(data, club_cutoff, Mgp_cutoff){
  mtx0 = data@images$slice1@coordinates
  mtx0 <- cbind(mtx0, t(data@assays$RCTD@data))
  mtx0 <- cbind(mtx0, data@meta.data)
  mtx0$Mgp <- data@assays$SCT@data['Mgp',]
  airCells <-  data.matrix(mtx0[mtx0$`Club cells` > quantile(mtx0$`Club cells`, club_cutoff), c("imagerow", "imagecol")])
  Th17 <- data.matrix(mtx0[mtx0$Mgp < quantile(mtx0$Mgp, Mgp_cutoff), c("imagerow", "imagecol", 'Th17')])
  Th1 <- data.matrix(mtx0[mtx0$Mgp < quantile(mtx0$Mgp, Mgp_cutoff), c("imagerow", "imagecol", 'Th1')])
  Th17.dist <- computeDist(Th17, airCells)
  Th1.dist <- computeDist(Th1, airCells)
  Th17.dist <- sum(Th17.dist*Th17[,3])/sum(Th17[,3])
  Th1.dist <- sum(Th1.dist*Th1[,3])/sum(Th1[,3])
  return(c(Th17.dist, Th1.dist))
}

th17_th1_dist_cutoff <- function(data, club_cutoff, Mgp_cutoff, cutoff){
  mtx0 = data@images$slice1@coordinates
  mtx0 <- cbind(mtx0, t(data@assays$RCTD@data))
  mtx0 <- cbind(mtx0, data@meta.data)
  mtx0$Mgp <- data@assays$SCT@data['Mgp',]
  airCells <-  data.matrix(mtx0[mtx0$`Club cells` > quantile(mtx0$`Club cells`, club_cutoff), c("imagerow", "imagecol")])
  Th17 <- data.matrix(mtx0[mtx0$Mgp < quantile(mtx0$Mgp, Mgp_cutoff), c("imagerow", "imagecol", 'Th17')])
  Th1 <- data.matrix(mtx0[mtx0$Mgp < quantile(mtx0$Mgp, Mgp_cutoff), c("imagerow", "imagecol", 'Th1')])
  Th17.dist <- computeDist(Th17, airCells)
  Th1.dist <- computeDist(Th1, airCells)
  Th17.dist <- Th17.dist[Th17.dist < cutoff]
  Th1.dist <- Th1.dist[Th1.dist < cutoff]
  Th17 <- Th17[names(Th17.dist),]
  Th1 <- Th1[names(Th1.dist),]
  Th17.dist <- sum(Th17.dist*Th17[,3])/sum(Th17[,3])
  Th1.dist <- sum(Th1.dist*Th1[,3])/sum(Th1[,3])
  return(c(Th17.dist, Th1.dist))
}

load('../Spatial/Rfiles/lung_spatial.RData')

th17_th1_dist(data = lung1, club_cutoff = 0.95, Mgp_cutoff = 0.95)
th17_th1_dist(data = lung2, club_cutoff = 0.95, Mgp_cutoff = 0.95)
th17_th1_dist(data = lung3, club_cutoff = 0.95, Mgp_cutoff = 0.95)
th17_th1_dist(data = lung4, club_cutoff = 0.95, Mgp_cutoff = 0.95)

lung1_dist <- data.frame('Th17' = NA, 'Th1' = NA, club_cutoff = rep(seq(0.90, 0.95, 0.01), each = 6), Mgp_cutoff = rep(seq(0.90, 0.95, 0.01), times = 6))
lung1_dist[,1:2] <- matrix(unlist(mclapply(1:36, function(i){th17_th1_dist(data = lung1, club_cutoff = lung1_dist[i,3], Mgp_cutoff = lung1_dist[i,4])}, mc.cores = 4)), ncol = 2, byrow = T)

lung2_dist <- data.frame('Th17' = NA, 'Th1' = NA, club_cutoff = rep(seq(0.90, 0.95, 0.01), each = 6), Mgp_cutoff = rep(seq(0.90, 0.95, 0.01), times = 6))
lung2_dist[,1:2] <- matrix(unlist(mclapply(1:36, function(i){th17_th1_dist(data = lung2, club_cutoff = lung2_dist[i,3], Mgp_cutoff = lung2_dist[i,4])}, mc.cores = 4)), ncol = 2, byrow = T)

lung3_dist <- data.frame('Th17' = NA, 'Th1' = NA, club_cutoff = rep(seq(0.90, 0.95, 0.01), each = 6), Mgp_cutoff = rep(seq(0.90, 0.95, 0.01), times = 6))
lung3_dist[,1:2] <- matrix(unlist(mclapply(1:36, function(i){th17_th1_dist(data = lung3, club_cutoff = lung3_dist[i,3], Mgp_cutoff = lung3_dist[i,4])}, mc.cores = 4)), ncol = 2, byrow = T)

lung4_dist <- data.frame('Th17' = NA, 'Th1' = NA, club_cutoff = rep(seq(0.90, 0.95, 0.01), each = 6), Mgp_cutoff = rep(seq(0.90, 0.95, 0.01), times = 6))
lung4_dist[,1:2] <- matrix(unlist(mclapply(1:36, function(i){th17_th1_dist(data = lung4, club_cutoff = lung4_dist[i,3], Mgp_cutoff = lung4_dist[i,4])}, mc.cores = 4)), ncol = 2, byrow = T)

save(lung1_dist, lung2_dist, lung3_dist, lung4_dist, file = 'distance.RData')

library(reshape2)

pdf('distance_A1.pdf', width = 4, height = 4)
tmp <- melt(lung1_dist, id.vars = c('club_cutoff', 'Mgp_cutoff'))
colnames(tmp)[3:4] <- c('Type', 'Value')
tmp$club_cutoff <- as.factor(tmp$club_cutoff)
tmp$Mgp_cutoff <- as.factor(tmp$Mgp_cutoff)
p <- ggpaired(tmp, x = "Type", y = "Value", color = "Type", line.color = "gray", line.size = 0.4)+
  rremove("x.text")
facet(p, facet.by = 'club_cutoff', ncol = 6)+
  stat_compare_means(paired = TRUE, method = 't.test', label = 'p.signif')+
  xlab('Cut-off of proportion of club cells (quantile)')+
  ylab('Weighted distance to airway (micrometer)')+
  ggtitle('Slice A1')+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf('distance_A2.pdf', width = 4, height = 4)
tmp <- melt(lung2_dist, id.vars = c('club_cutoff', 'Mgp_cutoff'))
colnames(tmp)[3:4] <- c('Type', 'Value')
tmp$club_cutoff <- as.factor(tmp$club_cutoff)
tmp$Mgp_cutoff <- as.factor(tmp$Mgp_cutoff)
p <- ggpaired(tmp, x = "Type", y = "Value", color = "Type", line.color = "gray", line.size = 0.4)+
  rremove("x.text")
facet(p, facet.by = 'club_cutoff', ncol = 6)+
  stat_compare_means(paired = TRUE, method = 't.test', label = 'p.signif')+
  xlab('Cut-off of proportion of club cells (quantile)')+
  ylab('Weighted distance to airway (micrometer)')+
  ggtitle('Slice A2')+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf('distance_A3.pdf', width = 4, height = 4)
tmp <- melt(lung3_dist, id.vars = c('club_cutoff', 'Mgp_cutoff'))
colnames(tmp)[3:4] <- c('Type', 'Value')
tmp$club_cutoff <- as.factor(tmp$club_cutoff)
tmp$Mgp_cutoff <- as.factor(tmp$Mgp_cutoff)
p <- ggpaired(tmp, x = "Type", y = "Value", color = "Type", line.color = "gray", line.size = 0.4)+
  rremove("x.text")
facet(p, facet.by = 'club_cutoff', ncol = 6)+
  stat_compare_means(paired = TRUE, method = 't.test', label = 'p.signif')+
  xlab('Cut-off of proportion of club cells (quantile)')+
  ylab('Weighted distance to airway (micrometer)')+
  ggtitle('Slice A3')+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf('distance_A4.pdf', width = 4, height = 4)
tmp <- melt(lung4_dist, id.vars = c('club_cutoff', 'Mgp_cutoff'))
colnames(tmp)[3:4] <- c('Type', 'Value')
tmp$club_cutoff <- as.factor(tmp$club_cutoff)
tmp$Mgp_cutoff <- as.factor(tmp$Mgp_cutoff)
p <- ggpaired(tmp, x = "Type", y = "Value", color = "Type", line.color = "gray", line.size = 0.4)+
  rremove("x.text")
facet(p, facet.by = 'club_cutoff', ncol = 6)+
  stat_compare_means(paired = TRUE, method = 't.test', label = 'p.signif')+
  xlab('Cut-off of proportion of club cells (quantile)')+
  ylab('Weighted distance to airway (micrometer)')+
  ggtitle('Slice A4')+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

lung1_dist_cutoff <- data.frame('Th17' = NA, 'Th1' = NA, club_cutoff = rep(seq(0.90, 0.95, 0.01), each = 6), Mgp_cutoff = rep(seq(0.90, 0.95, 0.01), times = 6))
lung1_dist_cutoff[,1:2] <- matrix(unlist(mclapply(1:36, function(i){th17_th1_dist_cutoff(data = lung1, club_cutoff = lung1_dist_cutoff[i,3], Mgp_cutoff = lung1_dist_cutoff[i,4], 1000)}, mc.cores = 4)), ncol = 2, byrow = T)

lung2_dist_cutoff <- data.frame('Th17' = NA, 'Th1' = NA, club_cutoff = rep(seq(0.90, 0.95, 0.01), each = 6), Mgp_cutoff = rep(seq(0.90, 0.95, 0.01), times = 6))
lung2_dist_cutoff[,1:2] <- matrix(unlist(mclapply(1:36, function(i){th17_th1_dist_cutoff(data = lung2, club_cutoff = lung2_dist_cutoff[i,3], Mgp_cutoff = lung2_dist_cutoff[i,4], 1000)}, mc.cores = 4)), ncol = 2, byrow = T)

lung3_dist_cutoff <- data.frame('Th17' = NA, 'Th1' = NA, club_cutoff = rep(seq(0.90, 0.95, 0.01), each = 6), Mgp_cutoff = rep(seq(0.90, 0.95, 0.01), times = 6))
lung3_dist_cutoff[,1:2] <- matrix(unlist(mclapply(1:36, function(i){th17_th1_dist_cutoff(data = lung3, club_cutoff = lung3_dist_cutoff[i,3], Mgp_cutoff = lung3_dist_cutoff[i,4], 1000)}, mc.cores = 4)), ncol = 2, byrow = T)

lung4_dist_cutoff <- data.frame('Th17' = NA, 'Th1' = NA, club_cutoff = rep(seq(0.90, 0.95, 0.01), each = 6), Mgp_cutoff = rep(seq(0.90, 0.95, 0.01), times = 6))
lung4_dist_cutoff[,1:2] <- matrix(unlist(mclapply(1:36, function(i){th17_th1_dist_cutoff(data = lung4, club_cutoff = lung4_dist_cutoff[i,3], Mgp_cutoff = lung4_dist_cutoff[i,4], 1000)}, mc.cores = 4)), ncol = 2, byrow = T)

save(lung1_dist_cutoff, lung2_dist_cutoff, lung3_dist_cutoff, lung4_dist_cutoff, file = 'distance_cutoff.RData')

pdf('distance_A1_cutoff.pdf', width = 4, height = 4)
tmp <- melt(lung1_dist_cutoff, id.vars = c('club_cutoff', 'Mgp_cutoff'))
colnames(tmp)[3:4] <- c('Type', 'Value')
tmp$club_cutoff <- as.factor(tmp$club_cutoff)
tmp$Mgp_cutoff <- as.factor(tmp$Mgp_cutoff)
p <- ggpaired(tmp, x = "Type", y = "Value", color = "Type", line.color = "gray", line.size = 0.4)+
  rremove("x.text")
facet(p, facet.by = 'club_cutoff', ncol = 6)+
  stat_compare_means(paired = TRUE, method = 't.test', label = 'p.signif')+
  xlab('Cut-off of proportion of club cells (quantile)')+
  ylab('Weighted distance to airway (micrometer)')+
  ggtitle('Slice A1')+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf('distance_A2_cutoff.pdf', width = 4, height = 4)
tmp <- melt(lung2_dist_cutoff, id.vars = c('club_cutoff', 'Mgp_cutoff'))
colnames(tmp)[3:4] <- c('Type', 'Value')
tmp$club_cutoff <- as.factor(tmp$club_cutoff)
tmp$Mgp_cutoff <- as.factor(tmp$Mgp_cutoff)
p <- ggpaired(tmp, x = "Type", y = "Value", color = "Type", line.color = "gray", line.size = 0.4)+
  rremove("x.text")
facet(p, facet.by = 'club_cutoff', ncol = 6)+
  stat_compare_means(paired = TRUE, method = 't.test', label = 'p.signif')+
  xlab('Cut-off of proportion of club cells (quantile)')+
  ylab('Weighted distance to airway (micrometer)')+
  ggtitle('Slice A2')+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf('distance_A3_cutoff.pdf', width = 4, height = 4)
tmp <- melt(lung3_dist_cutoff, id.vars = c('club_cutoff', 'Mgp_cutoff'))
colnames(tmp)[3:4] <- c('Type', 'Value')
tmp$club_cutoff <- as.factor(tmp$club_cutoff)
tmp$Mgp_cutoff <- as.factor(tmp$Mgp_cutoff)
p <- ggpaired(tmp, x = "Type", y = "Value", color = "Type", line.color = "gray", line.size = 0.4)+
  rremove("x.text")
facet(p, facet.by = 'club_cutoff', ncol = 6)+
  stat_compare_means(paired = TRUE, method = 't.test', label = 'p.signif')+
  xlab('Cut-off of proportion of club cells (quantile)')+
  ylab('Weighted distance to airway (micrometer)')+
  ggtitle('Slice A3')+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf('distance_A4_cutoff.pdf', width = 4, height = 4)
tmp <- melt(lung4_dist_cutoff, id.vars = c('club_cutoff', 'Mgp_cutoff'))
colnames(tmp)[3:4] <- c('Type', 'Value')
tmp$club_cutoff <- as.factor(tmp$club_cutoff)
tmp$Mgp_cutoff <- as.factor(tmp$Mgp_cutoff)
p <- ggpaired(tmp, x = "Type", y = "Value", color = "Type", line.color = "gray", line.size = 0.4)+
  rremove("x.text")
facet(p, facet.by = 'club_cutoff', ncol = 6)+
  stat_compare_means(paired = TRUE, method = 't.test', label = 'p.signif')+
  xlab('Cut-off of proportion of club cells (quantile)')+
  ylab('Weighted distance to airway (micrometer)')+
  ggtitle('Slice A4')+
  theme(plot.title = element_text(hjust = 0.5))
dev.off()
#######

th17_th1_dist_true <- function(data){
  mtx0 = data@images$slice1@coordinates
  mtx0 <- cbind(mtx0, t(data@assays$RCTD@data))
  mtx0 <- cbind(mtx0, data@meta.data)
  mtx0$Mgp <- data@assays$SCT@data['Mgp',]
  airCells <-  data.matrix(mtx0[mtx0$airway == 1, c("imagerow", "imagecol")])
  Th17 <- data.matrix(mtx0[mtx0$vessel == 0, c("imagerow", "imagecol", 'Th17')])
  Th1 <- data.matrix(mtx0[mtx0$vessel == 0, c("imagerow", "imagecol", 'Th1')])
  Th17.dist <- computeDist(Th17, airCells)
  Th1.dist <- computeDist(Th1, airCells)
  Th17.dist <- sum(Th17.dist*Th17[,3])/sum(Th17[,3])
  Th1.dist <- sum(Th1.dist*Th1[,3])/sum(Th1[,3])
  return(c(Th17.dist, Th1.dist))
}

th17_th1_dist_true_cutoff <- function(data, cutoff){
  mtx0 = data@images$slice1@coordinates
  mtx0 <- cbind(mtx0, t(data@assays$RCTD@data))
  mtx0 <- cbind(mtx0, data@meta.data)
  mtx0$Mgp <- data@assays$SCT@data['Mgp',]
  airCells <-  data.matrix(mtx0[mtx0$airway == 1, c("imagerow", "imagecol")])
  Th17 <- data.matrix(mtx0[mtx0$vessel == 0, c("imagerow", "imagecol", 'Th17')])
  Th1 <- data.matrix(mtx0[mtx0$vessel == 0, c("imagerow", "imagecol", 'Th1')])
  Th17.dist <- computeDist(Th17, airCells)
  Th1.dist <- computeDist(Th1, airCells)
  Th17.dist <- Th17.dist[Th17.dist < cutoff]
  Th1.dist <- Th1.dist[Th1.dist < cutoff]
  Th17 <- Th17[names(Th17.dist),]
  Th1 <- Th1[names(Th1.dist),]
  Th17.dist <- sum(Th17.dist*Th17[,3])/sum(Th17[,3])
  Th1.dist <- sum(Th1.dist*Th1[,3])/sum(Th1[,3])
  return(c(Th17.dist, Th1.dist))
}

check_contour <- function(type){
  p1 <- SpatialFeaturePlot(lung1, type, pt.size.factor = 1, alpha = c(0,1))
  p2 <- SpatialFeaturePlot(lung2, type, pt.size.factor = 1, alpha = c(0,1))
  p3 <- SpatialFeaturePlot(lung3, type, pt.size.factor = 1, alpha = c(0,1))
  p4 <- SpatialFeaturePlot(lung4, type, pt.size.factor = 1, alpha = c(0,1))
  plot <- (p1 | p2 | p3 | p4)
  return(plot)
}

lung_dist <- data.frame('Th17' = NA, 'Th1' = NA, 'Slice' = paste0('A', 1:4), 'Condition' = rep(c('Immunized Mouse', 'Re-challenged Mouse'), each = 2))
lung_dist[1,1:2] <- th17_th1_dist_true(data = lung1)
lung_dist[2,1:2] <- th17_th1_dist_true(data = lung2)
lung_dist[3,1:2] <- th17_th1_dist_true(data = lung3)
lung_dist[4,1:2] <- th17_th1_dist_true(data = lung4)

pdf('distance_true.pdf', width = 4, height = 4)
tmp <- melt(lung_dist, id.vars = c('Slice', 'Condition'))
colnames(tmp)[3:4] <- c('Type', 'Value')
p <- ggpaired(tmp, x = "Type", y = "Value", color = "Slice", line.color = "gray", line.size = 0.4)
facet(p, facet.by = 'Condition', ncol = 4)+
  xlab('Conditions')+
  ylab('Weighted distance to airway (micrometer)')
dev.off()

lung_dist <- data.frame('Th17' = NA, 'Th1' = NA, 'Slice' = paste0('A', 1:4), 'Condition' = rep(c('Immunized Mouse', 'Re-challenged Mouse'), each = 2))
lung_dist[1,1:2] <- th17_th1_dist_true_cutoff(data = lung1, 1000)
lung_dist[2,1:2] <- th17_th1_dist_true_cutoff(data = lung2, 1000)
lung_dist[3,1:2] <- th17_th1_dist_true_cutoff(data = lung3, 1000)
lung_dist[4,1:2] <- th17_th1_dist_true_cutoff(data = lung4, 1000)

pdf('distance_true_cutoff.pdf', width = 4, height = 4)
tmp <- melt(lung_dist, id.vars = c('Slice', 'Condition'))
colnames(tmp)[3:4] <- c('Type', 'Value')
p <- ggpaired(tmp, x = "Type", y = "Value", color = "Slice", line.color = "gray", line.size = 0.4)
facet(p, facet.by = 'Condition', ncol = 4)+
  xlab('Conditions')+
  ylab('Weighted distance to airway (micrometer)')
dev.off()

pdf('contour.pdf', width = 32, height = 16)
check_contour('airway') / check_contour('vessel')
dev.off()

th17_th1_dist_true_cutoff_include <- function(data, cutoff){
  mtx0 = data@images$slice1@coordinates
  mtx0 <- cbind(mtx0, t(data@assays$RCTD@data))
  mtx0 <- cbind(mtx0, data@meta.data)
  mtx0$Mgp <- data@assays$SCT@data['Mgp',]
  airCells <-  data.matrix(mtx0[mtx0$airway == 1, c("imagerow", "imagecol")])
  Th17 <- data.matrix(mtx0[, c("imagerow", "imagecol", 'Th17')])
  Th1 <- data.matrix(mtx0[, c("imagerow", "imagecol", 'Th1')])
  Th17.dist <- computeDist(Th17, airCells)
  Th1.dist <- computeDist(Th1, airCells)
  Th17.dist <- Th17.dist[Th17.dist < cutoff]
  Th1.dist <- Th1.dist[Th1.dist < cutoff]
  Th17 <- Th17[names(Th17.dist),]
  Th1 <- Th1[names(Th1.dist),]
  Th17.dist <- sum(Th17.dist*Th17[,3])/sum(Th17[,3])
  Th1.dist <- sum(Th1.dist*Th1[,3])/sum(Th1[,3])
  return(c(Th17.dist, Th1.dist))
}

lung_dist <- data.frame('Th17' = NA, 'Th1' = NA, 'Slice' = paste0('A', 1:4), 'Condition' = rep(c('Immunized Mouse', 'Re-challenged Mouse'), each = 2))
lung_dist[1,1:2] <- th17_th1_dist_true_cutoff_include(data = lung1, 1000)
lung_dist[2,1:2] <- th17_th1_dist_true_cutoff_include(data = lung2, 1000)
lung_dist[3,1:2] <- th17_th1_dist_true_cutoff_include(data = lung3, 1000)
lung_dist[4,1:2] <- th17_th1_dist_true_cutoff_include(data = lung4, 1000)

pdf('distance_true_cutoff_include.pdf', width = 4, height = 4)
tmp <- melt(lung_dist, id.vars = c('Slice', 'Condition'))
colnames(tmp)[3:4] <- c('Type', 'Value')
p <- ggpaired(tmp, x = "Type", y = "Value", color = "Slice", line.color = "gray", line.size = 0.4)
facet(p, facet.by = 'Condition', ncol = 4)+
  xlab('Conditions')+
  ylab('Weighted distance to airway (micrometer)')
dev.off()