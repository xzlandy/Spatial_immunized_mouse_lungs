library(Seurat)
library(tiff)

setwd('~/Box/RWorkSpace/Spatial_combine/Spatial/Rfiles/')

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
  return(data)
}

contour_vessel <- function(data, tif, cutoff){
  vessel <- tif[[2]]
  dim(vessel)
  rownames(vessel) <- 1:nrow(vessel)
  colnames(vessel) <- 1:ncol(vessel)
  
  # data@images$vessel <- data@images$slice1
  # data@images$vessel@image <- vessel
  
  mtx0 = data@images$slice1@coordinates
  scale.factor <- data@images$slice1@scale.factors$lowres
  mtx0$tifrow <- round(mtx0$imagerow * scale.factor)
  mtx0$tifcol <- round(mtx0$imagecol * scale.factor)
  mtx0$vessel <- sapply(1:nrow(mtx0), function(i){vessel[mtx0$tifrow[i], mtx0$tifcol[i]]})
  data$vessel <- ifelse(mtx0$vessel > cutoff, 1, 0)
  print(SpatialFeaturePlot(data, 'vessel', images = 'slice1', pt.size.factor = 1, alpha = c(0,1)))
  return(data)
}

contour_airway <- function(data, cutoff){
  data0 <- as.data.frame(t(data.matrix(data@assays$RCTD@data)))
  data$airway <- ifelse(data0$`Club cells` > cutoff, 1, 0)
  print(SpatialFeaturePlot(data, 'airway', images = 'slice1', pt.size.factor = 1, alpha = c(0,1)))
  return(data)
}

lung1 <- read_spatial_RCTD(spatial = 'A1_KPx2_12min.rds', rctd = '../../Deconvolution/spatial_a1_RCTD_transferred_rna_reference.RData')
lung2 <- read_spatial_RCTD(spatial = 'A2_KPx2_6min.rds', rctd = '../../Deconvolution/spatial_a2_RCTD_transferred_rna_reference.RData')
lung3 <- read_spatial_RCTD(spatial = 'A3_KPx3_12min.rds', rctd = '../../Deconvolution/spatial_a3_RCTD_transferred_rna_reference.RData')
lung4 <- read_spatial_RCTD(spatial = 'A4_KPx3_6min.rds', rctd = '../../Deconvolution/spatial_a4_RCTD_transferred_rna_reference.RData')

tif1 <- readTIFF('../../QuPath/A1_blood_vessel.ome_resized.tif', all = T)
tif2 <- readTIFF('../../QuPath/A2_blood_vessel.ome_resized.tif', all = T)
tif3 <- readTIFF('../../QuPath/A3_blood_vessel.ome_resized.tif', all = T)
tif4 <- readTIFF('../../QuPath/A4_blood_vessel.ome_resized.tif', all = T)

lung1 <- contour_vessel(lung1, tif1, 0.50)
lung2 <- contour_vessel(lung2, tif2, 0.50)
lung3 <- contour_vessel(lung3, tif3, 0.50)
lung4 <- contour_vessel(lung4, tif4, 0.50)

lung1 <- contour_airway(lung1, 0.20)
lung2 <- contour_airway(lung2, 0.20)
lung3 <- contour_airway(lung3, 0.1)
lung4 <- contour_airway(lung4, 0.1)

save(lung1, lung2, lung3, lung4, file = 'lung_spatial.RData')

tmp <- merge(lung1, lung2, add.cell.ids = c('A1', 'A2'))
tmp1 <- merge(lung3, lung4, add.cell.ids = c('A3', 'A4'))
merge.lung <- merge(tmp, tmp1)
rm(tmp, tmp1)

names(merge.lung@images) <- c("slice_A1", "slice_A2", "slice_A3", "slice_A4")

DefaultAssay(merge.lung) <- "SCT"
VariableFeatures(merge.lung) <- c(VariableFeatures(lung1, assay = 'SCT'), VariableFeatures(lung2, assay = 'SCT'), VariableFeatures(lung3, assay = 'SCT'), VariableFeatures(lung4, assay = 'SCT'))
merge.lung <- RunPCA(merge.lung, verbose = FALSE)
merge.lung <- FindNeighbors(merge.lung, dims = 1:30)
merge.lung <- FindClusters(merge.lung, verbose = FALSE)
merge.lung <- RunUMAP(merge.lung, dims = 1:30)

save(merge.lung, file = 'merge_lung_spatial.RData')

DefaultAssay(merge.lung) <- "RCTD"

pdf('spatial_predicted_celltypes_1.pdf', width = 16, height = 12)
SpatialFeaturePlot(merge.lung, rownames(merge.lung)[1:3])
dev.off()

pdf('spatial_predicted_celltypes_2.pdf', width = 16, height = 12)
SpatialFeaturePlot(merge.lung, rownames(merge.lung)[4:6])
dev.off()

pdf('spatial_predicted_celltypes_3.pdf', width = 16, height = 12)
SpatialFeaturePlot(merge.lung, rownames(merge.lung)[7:9])
dev.off()

pdf('spatial_predicted_celltypes_4.pdf', width = 16, height = 12)
SpatialFeaturePlot(merge.lung, rownames(merge.lung)[10:12])
dev.off()

pdf('spatial_predicted_celltypes_5.pdf', width = 16, height = 12)
SpatialFeaturePlot(merge.lung, rownames(merge.lung)[13:15])
dev.off()

pdf('spatial_predicted_celltypes_all.pdf', width = 16, height = 60)
SpatialFeaturePlot(merge.lung, rownames(merge.lung), ncol = 4)
dev.off()
