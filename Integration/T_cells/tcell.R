library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)

setwd('~/Box/RWorkSpace/Spatial_combine/Integration/T_cells/')

load('../../Annotation/rna_reference.RData')
rna_reference <- subset(rna_reference, cell_type == 'Low quality cells', invert = T)
rna_reference$cell_type <- factor(rna_reference$cell_type, levels = c('Alveolar epithelial cells', 'Club cells', 'Fibroblast', 'Endothelial cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'T cells', 'NK cells', 'Erythrocytes'))

load('../annotated_atac.RData')
DefaultAssay(data) <- 'SCT'

p1 <- FeaturePlot(rna_reference, c('Il17a', 'Il17f', 'Il21', 'Il22', 'Rorc', 'Rora', 'Stat3'), ncol = 7)
p2 <- FeaturePlot(data, c('Il17a', 'Il17f', 'Il21', 'Il22', 'Rorc', 'Rora', 'Stat3'), ncol = 7)
pdf('Th17_markers.pdf', width = 28, height = 8)
p1 / p2
dev.off()

p1 <- FeaturePlot(rna_reference, c('Ifng', 'Tnf', 'Lta', 'Il2', 'Tbx21', 'Stat1', 'Stat4'), ncol = 7)
p2 <- FeaturePlot(data, c('Ifng', 'Tnf', 'Lta', 'Il2', 'Tbx21', 'Stat1', 'Stat4'), ncol = 7)
pdf('Th1_markers.pdf', width = 28, height = 8)
p1 / p2
dev.off()

tcell <- subset(data, predicted.id == 'T cells')
DefaultAssay(tcell) <- 'peaks'
tcell <- FindTopFeatures(tcell, min.cutoff = 'q0')
tcell <- RunTFIDF(tcell)
tcell <- RunSVD(tcell)

DepthCor(tcell)

tcell <- RunUMAP(object = tcell, reduction = 'lsi', dims = 2:40)
tcell <- FindNeighbors(object = tcell, reduction = 'lsi', dims = 2:40)
tcell <- FindClusters(object = tcell, algorithm = 3, resolution = 0.2, graph.name = "peaks_snn")

pdf('UMAP_clusters.pdf', width = 8, height = 8)
DimPlot(tcell, label = T) + NoLegend()
dev.off()

DefaultAssay(tcell) <- "ACTIVITY"
tcell <- SCTransform(tcell, assay = "ACTIVITY")

DefaultAssay(tcell) <- "SCT"
p1 <- FeaturePlot(tcell, c('Il17a', 'Il17f', 'Il21', 'Il22', 'Rorc', 'Rora', 'Stat3'), ncol = 7)
p2 <- FeaturePlot(tcell, c('Ifng', 'Tnf', 'Lta', 'Il2', 'Tbx21', 'Stat1', 'Stat4'), ncol = 7)
pdf('scatter_Th17_Th1.pdf', width = 28, height = 8)
p1 / p2
dev.off()

markers <- FindAllMarkers(tcell, logfc.threshold = 0.25, only.pos = T)
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

pdf('heatmap_clusters.pdf', width = 16, height = 16)
DoHeatmap(tcell, features = top20$gene) + NoLegend()
dev.off()

Idents(tcell) <- 'seurat_clusters'
DefaultAssay(tcell) <- 'SCT'
markers <- FindAllMarkers(tcell, logfc.threshold = 0.25, only.pos = T, test.use = 'MAST')
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

pdf('heatmap_clusters_mast.pdf', width = 16, height = 16)
DoHeatmap(tcell, features = top20$gene) + NoLegend()
dev.off()

DefaultAssay(tcell) <- 'peaks'
da_peaks <- FindMarkers(
  object = tcell,
  ident.1 = 0,
  ident.2 = 2,
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

open_Th17 <- rownames(da_peaks[da_peaks$avg_log2FC > 0.5, ])
open_Th1 <- rownames(da_peaks[da_peaks$avg_log2FC < -0.5, ])

closest_genes_Th17 <- ClosestFeature(tcell, regions = open_Th17)
closest_genes_Th1 <- ClosestFeature(tcell, regions = open_Th1)
closest_genes <- ClosestFeature(tcell, regions = rownames(da_peaks))

da_peaks <- cbind(da_peaks, closest_genes)

pdf('peaks_Ifng.pdf', width = 8, height = 4)
CoveragePlot(
  object = tcell,
  region = 'chr10-118445262-118445892',
  extend.upstream = 20000,
  extend.downstream = 160000
)
dev.off()

pdf('peaks_Rora.pdf', width = 8, height = 4)
CoveragePlot(
  object = tcell,
  region = 'chr9-69197590-69348131',
  extend.upstream = 20000,
  extend.downstream = 20000
)
dev.off()

pdf('peaks_Tbx21.pdf', width = 8, height = 4)
CoveragePlot(
  object = tcell,
  region = 'chr11-97102099-97114649',
  extend.upstream = 20000,
  extend.downstream = 20000
)
dev.off()

####################
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)

# add motif information
tcell <- AddMotifs(
  object = tcell,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

da_peaks <- FindMarkers(
  object = tcell,
  ident.1 = 0,
  ident.2 = 2,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)

da_peaks_r <- FindMarkers(
  object = tcell,
  ident.1 = 2,
  ident.2 = 0,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
top.da.peak_r <- rownames(da_peaks_r[da_peaks_r$p_val < 0.005, ])

# find peaks open in Pvalb or Sst cells
open.peaks <- AccessiblePeaks(tcell, idents = c(0,2))

# # match the overall GC content in the peak set
# meta.feature <- GetAssayData(tcell, assay = "peaks", slot = "meta.features")
# peaks.matched <- MatchRegionStats(
#   meta.feature = meta.feature[open.peaks, ],
#   query.feature = meta.feature[top.da.peak, ],
#   n = 10000
# )

# test enrichment
enriched.motifs <- FindMotifs(
  object = tcell,
  features = top.da.peak
)

enriched.motifs_r <- FindMotifs(
  object = tcell,
  features = top.da.peak_r
)

pdf('motif_Th17.pdf', width = 4, height = 4)
MotifPlot(
  object = tcell,
  motifs = head(rownames(enriched.motifs))
)
dev.off()

pdf('motif_Th1.pdf', width = 4, height = 4)
MotifPlot(
  object = tcell,
  motifs = head(rownames(enriched.motifs_r))
)
dev.off()

# gather the footprinting information for sets of motifs
tcell <- Footprint(
  object = tcell,
  motif.name = c("RORC", "RORA", "STAT3", "TBX21", "STAT1"),
  genome = BSgenome.Mmusculus.UCSC.mm10
)
# plot the footprint data for each group of cells
pdf('footprint.pdf', width = 8, height = 16)
p2 <- PlotFootprint(tcell, features = c("RORC", "RORA", "STAT3", "TBX21", "STAT1"))
p2 + patchwork::plot_layout(ncol = 1)
dev.off()

pdf('footprint_Th17_Th1.pdf', width = 8, height = 16)
p2 <- PlotFootprint(tcell, features = c("RORC", "RORA", "STAT3", "TBX21", "STAT1"), idents = c(0,2))
p2 + patchwork::plot_layout(ncol = 1)
dev.off()
############
tcell <- RunChromVAR(
  object = tcell,
  genome = BSgenome.Mmusculus.UCSC.mm10
)
DefaultAssay(tcell) <- 'chromvar'

pdf('motif_activity.pdf', width = 12, height = 8)
FeaturePlot(
  object = tcell,
  features = c("MA1151.1", "MA0071.1", 'MA0144.2', 'MA0690.1', 'MA0137.3'),
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  ncol = 3
)
dev.off()

differential.activity <- FindMarkers(
  object = tcell,
  ident.1 = 0,
  ident.2 = 2,
  only.pos = TRUE,
  test.use = 'LR',
  latent.vars = 'nCount_peaks'
)
MotifPlot(
  object = tcell,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)

tcell <- RenameIdents(tcell,
             '0' = 'Th17',
             '1' = 'Other T cells 1',
             '2' = 'Th1',
             '3' = 'Other T cells 2')
tcell$cell_type <- tcell@active.ident

data$cell_type <- as.character(data$predicted.id)
data@meta.data[rownames(tcell@meta.data[tcell$cell_type == 'Th17',]),]$cell_type <- 'Th17'
data@meta.data[rownames(tcell@meta.data[tcell$cell_type == 'Th1',]),]$cell_type <- 'Th1'
data@meta.data[rownames(tcell@meta.data[tcell$cell_type == 'Other T cells 1',]),]$cell_type <- 'Other T cells 1'
data@meta.data[rownames(tcell@meta.data[tcell$cell_type == 'Other T cells 2',]),]$cell_type <- 'Other T cells 2'
data$cell_type <- factor(data$cell_type, levels = c('Alveolar epithelial cells', 'Club cells', 'Fibroblast', 'Endothelial cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'Th17', 'Th1', 'Other T cells 1', 'Other T cells 2', 'NK cells', 'Erythrocytes'))

Idents(tcell) <- 'cell_type'
pdf('UMAP_cell_types_tcell.pdf', width = 6, height = 6)
DimPlot(tcell)
dev.off()

Idents(data) <- 'cell_type'
pdf('UMAP_cell_types.pdf', width = 8, height = 8)
DimPlot(data)
dev.off()

save(tcell, data, file = 'tcell.RData')

DefaultAssay(tcell) <- 'SCT'
markers <- FindAllMarkers(tcell, logfc.threshold = 0.25, only.pos = T)
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

pdf('heatmap_clusters_annotated.pdf', width = 10, height = 10)
DoHeatmap(tcell, features = top20$gene, label = F)
dev.off()

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
pdf('heatmap_clusters_annotated_downsample.pdf', width = 6, height = 6)
DoHeatmap(subset(tcell, downsample = 100), features = top10$gene, label = F)
dev.off()

Idents(tcell) <- 'cell_type'
DefaultAssay(tcell) <- 'peaks'
da_peaks <- FindMarkers(
  object = tcell,
  ident.1 = 'Th17',
  ident.2 = 'Th1',
  min.pct = 0.1,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

open_Th17 <- rownames(da_peaks[da_peaks$avg_log2FC > 0.5, ])
open_Th1 <- rownames(da_peaks[da_peaks$avg_log2FC < -0.5, ])

closest_genes_Th17 <- ClosestFeature(tcell, regions = open_Th17)
closest_genes_Th1 <- ClosestFeature(tcell, regions = open_Th1)
closest_genes <- ClosestFeature(tcell, regions = rownames(da_peaks))

da_peaks <- cbind(da_peaks, closest_genes)

pdf('peaks_Ifng_annotated.pdf', width = 6, height = 3)
CoveragePlot(
  object = tcell,
  region = 'Ifng',
  extend.upstream = 10000,
  extend.downstream = 160000
)
dev.off()

pdf('peaks_Rora_annotated.pdf', width = 6, height = 3)
CoveragePlot(
  object = tcell,
  region = 'Rora',
  extend.upstream = 10000,
  extend.downstream = 10000
)
dev.off()

pdf('peaks_Rorc_annotated.pdf', width = 6, height = 3)
CoveragePlot(
  object = tcell,
  region = 'Rorc',
  extend.upstream = 10000,
  extend.downstream = 10000
)
dev.off()

pdf('peaks_Tbx21_annotated.pdf', width = 6, height = 3)
CoveragePlot(
  object = tcell,
  region = 'Tbx21',
  extend.upstream = 10000,
  extend.downstream = 10000
)
dev.off()

pdf('peaks_Il17a_annotated.pdf', width = 6, height = 3)
CoveragePlot(
  object = tcell,
  region = 'Il17a',
  extend.upstream = 10000,
  extend.downstream = 10000
)
dev.off()

pdf('footprint_annotated.pdf', width = 12, height = 3)
p2 <- PlotFootprint(tcell, features = c("RORC", "TBX21"))
p2 + patchwork::plot_layout(ncol = 2)
dev.off()

pdf('footprint_Th17_Th1_annotated.pdf', width = 12, height = 3)
p2 <- PlotFootprint(tcell, features = c("RORC", "TBX21"), idents = c('Th17', 'Th1'))
p2 + patchwork::plot_layout(ncol = 2)
dev.off()


