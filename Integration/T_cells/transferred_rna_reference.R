library(Signac)
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)

setwd('~/Box/RWorkSpace/Spatial_combine/Integration/T_cells/')

load('../../Annotation/rna_reference.RData')
rna_reference <- subset(rna_reference, cell_type == 'Low quality cells', invert = T)
rna_reference$cell_type <- factor(rna_reference$cell_type, levels = c('Alveolar epithelial cells', 'Club cells', 'Fibroblast', 'Endothelial cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'T cells', 'NK cells', 'Erythrocytes'))

load('tcell.RData')

rna_tcell <- subset(rna_reference, cell_type == 'T cells')
rna_tcell <- SCTransform(rna_tcell, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
rna_tcell <- RunPCA(rna_tcell, verbose = FALSE)
ElbowPlot(rna_tcell, ndims = 50)

rna_tcell <- RunUMAP(rna_tcell, dims = 1:50, verbose = FALSE)

rna_tcell <- FindNeighbors(rna_tcell, dims = 1:50, verbose = FALSE)
rna_tcell <- FindClusters(rna_tcell, resolution = 1)

pdf('UMAP_clusters_rna.pdf', width = 8, height = 8)
DimPlot(rna_tcell, label = TRUE) + NoLegend()
dev.off()

transfer.anchors <- FindTransferAnchors(reference = tcell, query = rna_tcell, 
                                        reference.assay = "SCT", query.assay = "SCT", reduction = "cca", normalization.method = 'SCT')

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = tcell$cell_type, 
                                     weight.reduction = rna_tcell[['pca']], dims = 1:50)

rna_tcell <- AddMetaData(rna_tcell, metadata = celltype.predictions)
Idents(rna_tcell) <- 'predicted.id'
pdf('UMAP_predicted_cell_types_rna.pdf', width = 6, height = 6)
DimPlot(rna_tcell)
dev.off()

p1 <- FeaturePlot(rna_tcell, c('Il17a', 'Il17f', 'Il21', 'Il22', 'Rorc', 'Rora', 'Stat3', 'Cd4', 'Ccr6', 'Ccr4', 'Klrb1c'), ncol = 11)
p2 <- FeaturePlot(rna_tcell, c('Ifng', 'Tnf', 'Lta', 'Il2', 'Tbx21', 'Stat1', 'Stat4', 'Cd4', 'Cxcr3', 'Ccr5', 'Il12rb2'), ncol = 11)
pdf('scatter_Th17_Th1_rna.pdf', width = 44, height = 8)
p1 / p2
dev.off()

markers <- FindAllMarkers(rna_tcell, logfc.threshold = 0.25, only.pos = T)
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

pdf('heatmap_predicted_cell_types_rna.pdf', width = 10, height = 10)
DoHeatmap(rna_tcell, features = top20$gene, label = F)
dev.off()

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf('heatmap_predicted_cell_types_rna_downsample.pdf', width = 6, height = 6)
DoHeatmap(subset(rna_tcell,downsample = 100), features = top10$gene, label = F) + theme(axis.text.y = element_text(face = 'italic', hjust = 1))
dev.off()

pdf('dotplot_markers_rna.pdf', width = 4, height = 4)
p <- DotPlot(rna_tcell, features = c('Rorc', 'Rora', 'Ccr6', 'Ccr4',
                                     'Ifng', 'Tbx21', 'Cxcr3', 'Ccr5', 'Il12rb2'))
p + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face = 'italic', hjust = 1)) + xlab('') + ylab('')
dev.off()

rna_tcell$cell_type <- rna_tcell$predicted.id
rna_reference$cell_type <- as.character(rna_reference$cell_type)
rna_reference@meta.data[rownames(rna_tcell@meta.data[rna_tcell$cell_type == 'Th17',]),]$cell_type <- 'Th17'
rna_reference@meta.data[rownames(rna_tcell@meta.data[rna_tcell$cell_type == 'Th1',]),]$cell_type <- 'Th1'
rna_reference@meta.data[rownames(rna_tcell@meta.data[rna_tcell$cell_type == 'Other T cells 1',]),]$cell_type <- 'Other T cells 1'
rna_reference@meta.data[rownames(rna_tcell@meta.data[rna_tcell$cell_type == 'Other T cells 2',]),]$cell_type <- 'Other T cells 2'
rna_reference$cell_type <- factor(rna_reference$cell_type, levels = c('Alveolar epithelial cells', 'Club cells', 'Fibroblast', 'Endothelial cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'Th17', 'Th1', 'Other T cells 1', 'Other T cells 2', 'NK cells', 'Erythrocytes'))

Idents(rna_reference) <- 'cell_type'
pdf('UMAP_predicted_cell_types_rna_all.pdf', width = 6, height = 6)
DimPlot(rna_reference)
dev.off()

save(rna_reference, file = 'transferred_rna_reference.RData')


