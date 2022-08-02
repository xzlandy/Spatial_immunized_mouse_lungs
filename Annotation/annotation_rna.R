library(Seurat)
library(ggplot2)
library(sctransform)
library(dplyr)

setwd('~/Box/RWorkSpace/Spatial_combine/Annotation/')
### SCtransform (not mito) ###
# Load raw count data of scRNA-seq
rna <- Read10X_h5("../RNA/outs/filtered_feature_bc_matrix.h5")
dim(rna)
rna <- CreateSeuratObject(counts = rna, assay = 'RNA', project = 'KC_mouse')
dim(rna)

# store mitochondrial percentage in object meta data
rna <- PercentageFeatureSet(rna, pattern = "^mt-", col.name = "percent.mt")
# Visualize QC metrics as a violin plot
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# run sctransform
rna_reference <- SCTransform(rna, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
rna_reference <- RunPCA(rna_reference, verbose = FALSE)
ElbowPlot(rna_reference, ndims = 50)

rna_reference <- RunUMAP(rna_reference, dims = 1:50, verbose = FALSE)

rna_reference <- FindNeighbors(rna_reference, dims = 1:50, verbose = FALSE)
rna_reference <- FindClusters(rna_reference, resolution = 0.2)
Idents(rna_reference) <- 'seurat_clusters'

pdf('UMAP_clusters.pdf', width = 4, height = 4)
DimPlot(rna_reference, label = TRUE) + NoLegend()
dev.off()

pdf('UMAP_mito.pdf', width = 4, height = 4)
FeaturePlot(rna_reference, 'percent.mt')
dev.off()

rna.markers <- FindAllMarkers(rna_reference, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(openxlsx)
write.xlsx(rna.markers, 'markers.xlsx', colNames = T, rowNames = F)
rna.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- rna.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf('heatmap_clusters.pdf', width = 16, height = 16)
DoHeatmap(rna_reference, features = top10$gene) + NoLegend()
dev.off()

top5 <- rna.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
pdf('heatmap_clusters_downsample.pdf', width = 6, height = 6)
DoHeatmap(subset(rna_reference,downsample =100), features = top5$gene) + NoLegend() + theme(axis.text.y = element_text(face = 'italic', hjust = 1))
dev.off()

pdf('scatter_T_cell.pdf', width = 12, height = 8)
FeaturePlot(rna_reference, c('Cd3d', 'Cd3e', 'Cd3g', 'Cd4', 'Cd8a', 'Cd8b1'), ncol = 3)
dev.off()

pdf('scatter_B_cell.pdf', width = 12, height = 4)
FeaturePlot(rna_reference, c('Ptprc', 'Cd19', 'Cd22'), ncol = 3)
dev.off()

pdf('scatter_DC.pdf', width = 24, height = 12)
FeaturePlot(rna_reference, c('Cd80', 'Cd86', 'Itgam', 'Itgax', 'Cd40', 'Ptprc', 'Cd83', 'Itgae', 'Clec9a', 'Cx3cr1', 'Ly75', 'Adgre1', 'Fcgr1', 'Cadm1', 'Sirpa', 'Xcr1', 'Cst3'), ncol = 6)
dev.off()

pdf('scatter_NK.pdf', width = 4, height = 4)
FeaturePlot(rna_reference, c('Ncr1'), ncol = 1)
dev.off()

pdf('scatter_Monocytes.pdf', width = 12, height = 8)
FeaturePlot(rna_reference, c('Itgam', 'Ly6c1', 'Csf1r', 'Adgre1', 'Cd14'), ncol = 3)
dev.off()

pdf('scatter_Macrophages.pdf', width = 32, height = 12)
FeaturePlot(rna_reference, c('Cd80', 'Cd86', 'Ccr5', 'Itgam', 'Itgax', 'Cd14', 'Fut4', 'Cd68', 'Adgre1', 'Fcgr1', 'Fcgr2b', 'Fcgr3', 'Lgals3', 'Itgal', 'Lamp2', 'Lilrb4a', 'Csf1r', 'Cd33', 'Tlr2', 'Tlr4', 'Mrc1' ,'Cd68', 'Marco'), ncol = 8)
dev.off()

pdf('scatter_Neutrophils.pdf', width = 12, height = 12)
FeaturePlot(rna_reference, c('Itgam', 'Ptprc', 'Ly6g', 'Gsr', 'Itgb2', 'Pglyrp1', 'S100a8'), ncol = 3)
dev.off()

pdf('scatter_Endo.pdf', width = 12, height = 4)
FeaturePlot(rna_reference, c('Mcam', 'Vcam1', 'Pecam1'), ncol = 3)
dev.off()

pdf('scatter_Epi.pdf', width = 12, height = 4)
FeaturePlot(rna_reference, c('Epcam', 'Sftpc', 'Scgb1a1'), ncol = 3)
dev.off()

pdf('scatter_Fibro.pdf', width = 4, height = 4)
FeaturePlot(rna_reference, c('Col3a1'), ncol = 1)
dev.off()

pdf('scatter_Erythrocytes.pdf', width = 12, height = 8)
FeaturePlot(rna_reference, c('Hba-a1', 'Hba-a2', 'Hbb-bt', 'Hbb-bs'), ncol = 3)
dev.off()

rna_reference$cell_type <- rna_reference$seurat_clusters
rna_reference$cell_type <- gsub(3, 'T cells', rna_reference$cell_type)
rna_reference$cell_type <- gsub(5, 'B cells', rna_reference$cell_type)
rna_reference$cell_type <- gsub(8, 'NK cells', rna_reference$cell_type)

rna_reference$cell_type <- gsub(1, 'Low quality cells', rna_reference$cell_type)
rna_reference$cell_type <- gsub(0, 'Monocytes', rna_reference$cell_type)
rna_reference$cell_type <- gsub(2, 'Macrophages', rna_reference$cell_type)
rna_reference$cell_type <- gsub(4, 'Dendritic cells', rna_reference$cell_type)
rna_reference$cell_type <- gsub(7, 'Neutrophils', rna_reference$cell_type)
rna_reference$cell_type <- gsub(9, 'Erythrocytes', rna_reference$cell_type)

rna_epi_endo <- subset(rna_reference, seurat_clusters == 6)
rna_epi_endo <- SCTransform(rna_epi_endo, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
rna_epi_endo <- RunPCA(rna_epi_endo, verbose = FALSE)
ElbowPlot(rna_epi_endo, ndims = 50)

rna_epi_endo <- RunUMAP(rna_epi_endo, dims = 1:50, verbose = FALSE)

rna_epi_endo <- FindNeighbors(rna_epi_endo, dims = 1:50, verbose = FALSE)
rna_epi_endo <- FindClusters(rna_epi_endo, resolution = 0.6)

pdf('UMAP_clusters_epi_endo.pdf', width = 4, height = 4)
DimPlot(rna_epi_endo, label = TRUE) + NoLegend()
dev.off()

rna_epi_endo.markers <- FindAllMarkers(rna_epi_endo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
rna_epi_endo.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
rna_epi_endo.markers.top10 <- rna_epi_endo.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf('heatmap_clusters_epi_endo.pdf', width = 8, height = 8)
DoHeatmap(rna_epi_endo, features = rna_epi_endo.markers.top10$gene) + NoLegend()
dev.off()

top5 <- rna_epi_endo.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf('heatmap_clusters_epi_endo_downsample.pdf', width = 3.5, height = 3.5)
DoHeatmap(subset(rna_epi_endo,downsample =100), features = top5$gene) + NoLegend() + theme(axis.text.y = element_text(face = 'italic', hjust = 1))
dev.off() 

pdf('scatter_Endo_epi_endo.pdf', width = 12, height = 8)
FeaturePlot(rna_epi_endo, c('Mcam', 'Vcam1', 'Pecam1', 'Cdh5'), ncol = 3)
dev.off()

pdf('scatter_Epi_epi_endo.pdf', width = 12, height = 4)
FeaturePlot(rna_epi_endo, c('Epcam', 'Sftpc', 'Scgb1a1'), ncol = 3)
dev.off()

pdf('scatter_Fibro_epi_endo.pdf', width = 8, height = 4)
FeaturePlot(rna_epi_endo, c('Col3a1', 'Col1a2'), ncol = 2)
dev.off()

rna_epi_endo$cell_type <- rna_epi_endo$seurat_clusters
rna_epi_endo$cell_type <- gsub(0, 'Endothelial cells', rna_epi_endo$cell_type)
rna_epi_endo$cell_type <- gsub(1, 'Fibroblast', rna_epi_endo$cell_type)
rna_epi_endo$cell_type <- gsub(2, 'Alveolar epithelial cells', rna_epi_endo$cell_type)
rna_epi_endo$cell_type <- gsub(3, 'Club cells', rna_epi_endo$cell_type)

Idents(rna_epi_endo) <- 'cell_type'
pdf('UMAP_clusters_epi_endo_annotated.pdf', width = 6, height = 4)
DimPlot(rna_epi_endo, label = F)
dev.off()

rna_reference@meta.data[rownames(rna_epi_endo@meta.data[rna_epi_endo$seurat_clusters == 0,]),]$cell_type <- 'Endothelial cells'
rna_reference@meta.data[rownames(rna_epi_endo@meta.data[rna_epi_endo$seurat_clusters == 1,]),]$cell_type <- 'Fibroblast'
rna_reference@meta.data[rownames(rna_epi_endo@meta.data[rna_epi_endo$seurat_clusters == 2,]),]$cell_type <- 'Alveolar epithelial cells'
rna_reference@meta.data[rownames(rna_epi_endo@meta.data[rna_epi_endo$seurat_clusters == 3,]),]$cell_type <- 'Club cells'

rna_reference$cell_type <- factor(rna_reference$cell_type, levels = c('Alveolar epithelial cells', 'Club cells', 'Fibroblast', 'Endothelial cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'T cells', 'NK cells', 'Erythrocytes', 'Low quality cells'))

Idents(rna_reference) <- 'cell_type'
pdf('UMAP_clusters_annotated.pdf', width = 8, height = 8)
DimPlot(rna_reference, label = F)
dev.off()

save(rna_reference, file = 'rna_reference.RData')

load('rna_reference.RData')
rna_reference <- subset(rna_reference, cell_type == 'Low quality cells', invert = T)
rna_reference$cell_type <- factor(rna_reference$cell_type, levels = c('Alveolar epithelial cells', 'Club cells', 'Fibroblast', 'Endothelial cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'T cells', 'NK cells', 'Erythrocytes'))

pdf('dotplot_markers.pdf', width = 8, height = 8)
DotPlot(rna_reference, features = c('Sftpc', 'Scgb1a1', 'Col3a1', 'Cdh5', 'Cd14', 'Cd68', 'Cst3', 'S100a8', 'Cd19', 'Cd3d', 'Ncr1', 'Hba-a1')) + RotatedAxis()
dev.off()

pdf('UMAP_clusters_annotated_without_mito.pdf', width = 6, height = 6)
DimPlot(rna_reference, label = F)
dev.off()

pdf('dotplot_markers_final.pdf', width = 5, height = 10)
p <- DotPlot(rna_reference, features = c('Sftpa1', 'Sftpb', 'Sftpc', 'Sftpd',
                                         'Scgb1a1', 'Muc5b', 'Scgb3a1', 'Scgb3a2',
                                         'Col3a1', 'Col1a2', 'Col1a1', 'Mfap4',
                                         'Cdh5', 'Mcam', 'Vcam1', 'Pecam1',
                                         'Cd14', 'Itgam',
                                         'Itgax', 'Cd68', 'Mrc1', 'Marco',
                                         'Aif1', 'H2-DMb1', 'H2-Eb1', 'H2-Aa',
                                         'Gsr', 'Pglyrp1', 'S100a8', 'Ly6g',
                                         'Ms4a1', 'Cd79a', 'Igkc', 'Cd19',
                                         'Cd3d', 'Cd3e', 'Cd3g',
                                         'Ncr1', 'Nkg7',
                                         'Hbb-bt', 'Hba-a2'))
p + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face = 'italic', hjust = 1), legend.position = 'bottom', legend.direction = 'vertical') + xlab('') + ylab('')
dev.off()

rna.markers <- FindAllMarkers(rna_reference, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(openxlsx)
write.xlsx(rna.markers, 'markers_final.xlsx', colNames = T, rowNames = F)
rna.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top5 <- rna.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
pdf('heatmap_clusters_final.pdf', width = 8, height = 8)
DoHeatmap(rna_reference, features = top5$gene, label = F)
dev.off()

top5 <- rna.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
pdf('heatmap_clusters_final_downsample.pdf', width = 6, height = 6)
DoHeatmap(subset(rna_reference,downsample =100), features = top5$gene, label = F) + theme(axis.text.y = element_text(face = 'italic', hjust = 1))
dev.off()
