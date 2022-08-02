library(Seurat)
library(openxlsx)
library(limma)
library(parallel)
library(reshape2)
library(splines)
library(ggplot2)
library(CellChat)
library(patchwork)
library(future)
library(ComplexHeatmap)

setwd('~/Box/RWorkSpace/Spatial_combine/CellChat/ST/')

load('../../Spatial/Rfiles/lung_spatial.RData')
load('../../Spatial/Rfiles/merge_lung_spatial.RData')

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
  
  data <- raw_data[,raw_data$cells == 1]
  data$cell_types <- NA
  data$cell_types[data$Alveolar_epithelial_cells == 1] <- 'Alveolar epithelial cells'
  data$cell_types[data$Club_cells == 1] <- 'Club cells'
  data$cell_types[data$Monocytes == 1] <- 'Monocytes'
  data$cell_types[data$Macrophages == 1] <- 'Macrophages'
  data$cell_types[data$Dendritic_cells == 1] <- 'Dendritic cells'
  data$cell_types[data$Neutrophils == 1] <- 'Neutrophils'
  data$cell_types[data$B_cells == 1] <- 'B cells'
  data$cell_types[data$Th17 == 1] <- 'Th17'
  data$cell_types[data$Th1 == 1] <- 'Th1'
  data$cell_types[data$Other_T_cells_1 == 1] <- 'Other T cells 1'
  data$cell_types[data$Other_T_cells_2 == 1] <- 'Other T cells 2'
  data$cell_types[data$NK_cells == 1] <- 'NK cells'
  data$cell_types <- factor(data$cell_types, levels = c('Alveolar epithelial cells', 'Club cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'Th17', 'Th1', 'Other T cells 1', 'Other T cells 2', 'NK cells'))
  
  Idents(data) <- 'cell_types'
  return(data)
}

run_cellchat <- function(raw_data){
  data <- define_cell(raw_data)
  
  data.input <- GetAssayData(data, assay = "SCT", slot = "data") # normalized data matrix
  labels <- Idents(data)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  
  CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  
  dplyr::glimpse(CellChatDB$interaction)
  # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 8) # do parallel
  #> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
  #> in future (>= 1.13.0) when running R from RStudio, because it is
  #> considered unstable. Because of this, plan("multicore") will fall
  #> back to plan("sequential"), and plan("multiprocess") will fall back to
  #> plan("multisession") - not plan("multicore") as in the past. For more details,
  #> how to control forked processing or not, and how to silence this warning in
  #> future R sessions, see ?future::supportsMulticore
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI network (optional)
  cellchat <- projectData(cellchat, PPI.mouse)
  
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  # cellchat <- filterCommunication(cellchat, min.cells = 5)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}

cellchat1 <- run_cellchat(lung1)
cellchat2 <- run_cellchat(lung2)
cellchat3 <- run_cellchat(lung3)
cellchat4 <- run_cellchat(lung4)

run_cellchat_merge <- function(raw_data1, raw_data2, ids){
  data1 <- define_cell(raw_data1)
  data2 <- define_cell(raw_data2)
  data <- merge(data1, data2, add.cell.ids = ids)
  data$cell_types <- factor(data$cell_types, levels = c('Alveolar epithelial cells', 'Club cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'Th17', 'Th1', 'Other T cells 1', 'Other T cells 2', 'NK cells'))
  Idents(data) <- 'cell_types'
  
  data.input <- GetAssayData(data, assay = "SCT", slot = "data") # normalized data matrix
  labels <- Idents(data)
  meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  
  CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  
  dplyr::glimpse(CellChatDB$interaction)
  # CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 8) # do parallel
  #> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
  #> in future (>= 1.13.0) when running R from RStudio, because it is
  #> considered unstable. Because of this, plan("multicore") will fall
  #> back to plan("sequential"), and plan("multiprocess") will fall back to
  #> plan("multisession") - not plan("multicore") as in the past. For more details,
  #> how to control forked processing or not, and how to silence this warning in
  #> future R sessions, see ?future::supportsMulticore
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI network (optional)
  cellchat <- projectData(cellchat, PPI.mouse)
  
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  # cellchat <- filterCommunication(cellchat, min.cells = 5)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  return(cellchat)
}

cellchat_KPx2 <- run_cellchat_merge(lung1, lung2, c('lung1', 'lung2'))
cellchat_KPx3 <- run_cellchat_merge(lung3, lung4, c('lung3', 'lung4'))

# save(cellchat1, cellchat2, cellchat3, cellchat4, cellchat_KPx2, cellchat_KPx3, file = 'ST_cellchat.RData')
load('../../Spatial/Rfiles/merge_lung_spatial.RData')
load('../../Integration/T_cells/transferred_rna_reference.RData')

check_gene_only <- function(gene){
  SpatialFeaturePlot(merge.lung, gene)
}

load('../scRNA/scRNA_cellchat.RData')

cellchat_scRNA <- cellchat

load('../scATAC/scATAC_cellchat.RData')

cellchat_scATAC <- cellchat

load('ST_cellchat.RData')

cellchat_KPx2 <- netAnalysis_computeCentrality(cellchat_KPx2, slot.name = "netP")
cellchat_KPx3 <- netAnalysis_computeCentrality(cellchat_KPx3, slot.name = "netP")

object.list_merge <- list(A1 = cellchat1, A2 = cellchat2, A3 = cellchat3, A4 = cellchat4)
cellchat_merge <- mergeCellChat(object.list_merge, add.names = names(object.list_merge), cell.prefix = T)
cellchat_merge

object.list_merge_KPx2 <- list(A1 = cellchat1, A2 = cellchat2)
cellchat_merge_KPx2 <- mergeCellChat(object.list_merge_KPx2, add.names = names(object.list_merge_KPx2), cell.prefix = T)
cellchat_merge_KPx2

object.list_merge_KPx3 <- list(A3 = cellchat3, A4 = cellchat4)
cellchat_merge_KPx3 <- mergeCellChat(object.list_merge_KPx3, add.names = names(object.list_merge_KPx3), cell.prefix = T)
cellchat_merge_KPx3

object.list <- list(`Immunized Mouse` = cellchat_KPx2, `Re-challenged Mouse` = cellchat_KPx3)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)
cellchat

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1:2)) + theme(axis.text.x = element_text(angle = 90))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1:2), measure = "weight") + theme(axis.text.x = element_text(angle = 90))
pdf('bar_compare_count_weight.pdf', width = 4, height = 4)
gg1 + gg2
dev.off()

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
pdf('heatmap_compare_count_weight.pdf', width = 8, height = 4)
gg1 + gg2
dev.off()

gg1 <- netVisual_heatmap(cellchat, cluster.rows = T, cluster.cols = T)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight", cluster.rows = T, cluster.cols = T)
#> Do heatmap based on a merged object
pdf('heatmap_compare_count_cluster.pdf', width = 5, height = 5)
gg1
dev.off()

pdf('heatmap_compare_weight_cluster.pdf', width = 5, height = 5)
gg2
dev.off()

gg1 <- netVisual_heatmap(cellchat_KPx2, color.heatmap = 'Reds')
gg2 <- netVisual_heatmap(cellchat_KPx3, color.heatmap = 'Reds')
pdf('heatmap_count.pdf', width = 8, height = 4)
gg1 + gg2
dev.off()

pdf('circle_count.pdf', width = 16, height = 8)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, comparison = c(2,1)) + theme(legend.position = 'bottom', legend.direction = 'vertical')
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison = c(2,1)) + theme(legend.position = 'none')
pdf('rank.pdf', width = 6, height = 12)
gg1 + gg2
dev.off()

pdf('circle_CXCL.pdf', width = 16, height = 8)
pathways.show <- c("CXCL")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# Access all the signaling pathways showing significant communications
pathways.show.all <- union(cellchat@netP$`Immunized Mouse`$pathways, cellchat@netP$`Re-challenged Mouse`$pathways)

whether_include <- function(pathway){
  index <- 1:2
  decision <- c(pathway %in% cellchat@netP$`Immunized Mouse`$pathways, pathway %in% cellchat@netP$`Re-challenged Mouse`$pathways)
  index <- index[decision]
  return(index)
}

for (j in 1:length(pathways.show.all)) {
  num <- whether_include(pathways.show.all[j])
  num_l <- length(num)
  pdf(paste0('circle_', pathways.show.all[j], ".pdf"), width = 6*num_l, height = 6)
  pathways.show <- pathways.show.all[j]
  weight.max <- getMaxWeight(object.list[num], slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
  par(mfrow = c(1,num_l), xpd=TRUE)
  for (i in num) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
  }
  dev.off()
}

check_pathway <- function(cellchat, pathways.show){
  print(netAnalysis_contribution(cellchat, signaling = pathways.show))
  pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  for (i in 1:nrow(pairLR)) {
    netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = pairLR[i,], layout = "circle")
  }
}

check_pathway_print <- function(pathways.show){
  pdf(paste0('individual_A1_', pathways.show, '.pdf'))
  try(check_pathway(cellchat1, pathways.show))
  dev.off()
  
  pdf(paste0('individual_A2_', pathways.show, '.pdf'))
  try(check_pathway(cellchat2, pathways.show))
  dev.off()
  
  pdf(paste0('individual_A3_', pathways.show, '.pdf'))
  try(check_pathway(cellchat3, pathways.show))
  dev.off()
  
  pdf(paste0('individual_A4_', pathways.show, '.pdf'))
  try(check_pathway(cellchat4, pathways.show))
  dev.off()
  
  pdf(paste0('individual_KPx3_', pathways.show, '.pdf'))
  try(check_pathway(cellchat_KPx3, pathways.show))
  dev.off()
  
  pdf(paste0('individual_KPx2_', pathways.show, '.pdf'))
  try(check_pathway(cellchat_KPx2, pathways.show))
  dev.off()
  
  pdf(paste0('individual_scRNA_', pathways.show, '.pdf'))
  try(check_pathway(cellchat_scRNA, pathways.show))
  dev.off()
  
  pdf(paste0('individual_scATAC_', pathways.show, '.pdf'))
  try(check_pathway(cellchat_scATAC, pathways.show))
  dev.off()
}

check_pathway_print('IL1')

check_pathway_print('CXCL')

check_pathway_print('CCL')

check_pathway_print('IL6')

check_pathway_print('CSF3')

check_pathway_print('CSF')

check_pathway_print('TGFb')

check_pathway_print('IFN-II')

check_pathway_print('LT')

check_gene_only(c('Cxcl1', 'Cxcr2', 'rctd_Neutrophils'))

check_gene_only(c('Il6', 'Il6ra', 'Il6st'))

check_gene_only(c('Il22', 'rctd_Th17'))

cellchat_merge_KPx3@meta$datasets = factor(cellchat_merge_KPx3@meta$datasets, levels = c("A3", "A4")) # set factor level
pdf('exp_IL6_KPx3.pdf', width = 4, height = 4)
plotGeneExpression(cellchat_merge_KPx3, signaling = "IL6", split.by = "datasets", colors.ggplot = T)
dev.off()

pdf('exp_TGFb_KPx3.pdf', width = 4, height = 4)
plotGeneExpression(cellchat_merge_KPx3, signaling = "TGFb", split.by = "datasets", colors.ggplot = T)
dev.off()

pdf('exp_TGFb_KPx2.pdf', width = 4, height = 4)
plotGeneExpression(cellchat_merge_KPx2, signaling = "TGFb", split.by = "datasets", colors.ggplot = T)
dev.off()

pdf('exp_LT_KPx3.pdf', width = 4, height = 4)
plotGeneExpression(cellchat_merge_KPx3, signaling = "LT", split.by = "datasets", colors.ggplot = T)
dev.off()
