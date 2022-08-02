library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd('~/Box/RWorkSpace/Spatial_combine/Integration/')
### Pre-processing workflow
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "mm10"

counts <- Read10X_h5(filename = "../ATAC/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "../ATAC/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = '../ATAC/outs/fragments.tsv.gz',
  annotation = annotation
)

data <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks',
  project = 'Kong_mouse',
  meta.data = metadata
)

#Analysis of ATAC data

# compute nucleosome signal score per cell
data <- NucleosomeSignal(object = data)

# compute TSS enrichment score per cell
data <- TSSEnrichment(object = data, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100
data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments

data$high.tss <- ifelse(data$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(data, group.by = 'high.tss') + NoLegend()

data$nucleosome_group <- ifelse(data$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = data, group.by = 'nucleosome_group', region = 'chr1-1-10000000')

pdf('violin_QC.pdf', width = 40, height = 8)
VlnPlot(
  object = data,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
dev.off()

### Create a gene activity matrix
# compute gene activities
gene.activities <- GeneActivity(data)

# add the gene activity matrix to the Seurat object as a new assay
# add gene activities as a new assay
data[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

save(data, file = 'data.RData')

load('data.RData')
data <- subset(x = data,
       subset = pct_reads_in_peaks > 15 &
         blacklist_ratio < 0.05 &
         nucleosome_signal < 4 &
         TSS.enrichment > 2)

DefaultAssay(data) <- "peaks"
data <- FindTopFeatures(data, min.cutoff = 'q0')
data <- RunTFIDF(data)
data <- RunSVD(data)

DepthCor(data)

data <- RunUMAP(object = data, reduction = 'lsi', dims = 2:40)
data <- FindNeighbors(object = data, reduction = 'lsi', dims = 2:40)
data <- FindClusters(object = data, algorithm = 3, resolution = 0.2, graph.name = "peaks_snn")

pdf('UMAP_cluster.pdf', width = 8, height = 8)
DimPlot(object = data, label = TRUE) + NoLegend()
dev.off()

# normalize gene activities
DefaultAssay(data) <- "ACTIVITY"
data <- SCTransform(data, assay = "ACTIVITY")

# Identify anchors
load('../Annotation/rna_reference.RData')
rna_reference <- subset(rna_reference, cell_type == 'Low quality cells', invert = T)
rna_reference$cell_type <- factor(rna_reference$cell_type, levels = c('Alveolar epithelial cells', 'Club cells', 'Fibroblast', 'Endothelial cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'T cells', 'NK cells', 'Erythrocytes'))

transfer.anchors <- FindTransferAnchors(reference = rna_reference, query = data, 
                                        reference.assay = "SCT", query.assay = "SCT", reduction = "cca", normalization.method = 'SCT')

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rna_reference$cell_type, 
                                     weight.reduction = data[["lsi"]], dims = 2:40)

data <- AddMetaData(data, metadata = celltype.predictions)

plot1 <- DimPlot(
  object = rna_reference,
  label = F,
  repel = TRUE) + ggtitle('scRNA-seq')

data$predicted.id <- factor(data$predicted.id, levels = c('Alveolar epithelial cells', 'Club cells', 'Fibroblast', 'Endothelial cells', 'Monocytes', 'Macrophages', 'Dendritic cells', 'Neutrophils', 'B cells', 'T cells', 'NK cells', 'Erythrocytes'))

plot2 <- DimPlot(
  object = data,
  group.by = 'predicted.id',
  label = F,
  repel = TRUE) + ggtitle('scATAC-seq')

pdf('UMAP_predicted_cell_types.pdf', width = 16, height = 8)
plot1 + plot2
dev.off()

pdf('UMAP_cell_types.pdf', width = 8, height = 8)
DimPlot(
  object = data,
  group.by = 'predicted.id',
  label = F,
  repel = TRUE)
dev.off()

pdf('scatter_T_cell.pdf', width = 12, height = 8)
FeaturePlot(data, c('Cd3d', 'Cd3e', 'Cd3g', 'Cd4', 'Cd8a', 'Cd8b1'), ncol = 3)
dev.off()

pdf('scatter_B_cell.pdf', width = 12, height = 4)
FeaturePlot(data, c('Ptprc', 'Cd19', 'Cd22'), ncol = 3)
dev.off()

pdf('scatter_DC.pdf', width = 24, height = 12)
FeaturePlot(data, c('Cd80', 'Cd86', 'Itgam', 'Itgax', 'Cd40', 'Ptprc', 'Cd83', 'Itgae', 'Clec9a', 'Cx3cr1', 'Ly75', 'Adgre1', 'Fcgr1', 'Cadm1', 'Sirpa', 'Xcr1', 'Cst3'), ncol = 6)
dev.off()

pdf('scatter_NK.pdf', width = 4, height = 4)
FeaturePlot(data, c('Ncr1'), ncol = 1)
dev.off()

pdf('scatter_Monocytes.pdf', width = 12, height = 8)
FeaturePlot(data, c('Itgam', 'Ly6c1', 'Csf1r', 'Adgre1', 'Cd14'), ncol = 3)
dev.off()

pdf('scatter_Macrophages.pdf', width = 32, height = 12)
FeaturePlot(data, c('Cd80', 'Cd86', 'Ccr5', 'Itgam', 'Itgax', 'Cd14', 'Fut4', 'Cd68', 'Adgre1', 'Fcgr1', 'Fcgr2b', 'Fcgr3', 'Lgals3', 'Itgal', 'Lamp2', 'Lilrb4a', 'Csf1r', 'Cd33', 'Tlr2', 'Tlr4', 'Mrc1' ,'Cd68', 'Marco'), ncol = 8)
dev.off()

pdf('scatter_Neutrophils.pdf', width = 12, height = 12)
FeaturePlot(data, c('Itgam', 'Ptprc', 'Ly6g', 'Gsr', 'Itgb2', 'Pglyrp1', 'S100a8'), ncol = 3)
dev.off()

pdf('scatter_Endo.pdf', width = 12, height = 8)
FeaturePlot(data, c('Mcam', 'Vcam1', 'Pecam1', 'Cdh5'), ncol = 3)
dev.off()

pdf('scatter_Epi.pdf', width = 12, height = 4)
FeaturePlot(data, c('Epcam', 'Sftpc', 'Scgb1a1'), ncol = 3)
dev.off()

pdf('scatter_Fibro.pdf', width = 8, height = 4)
FeaturePlot(data, c('Col3a1', 'Col1a2'), ncol = 2)
dev.off()

pdf('scatter_Erythrocytes.pdf', width = 12, height = 8)
FeaturePlot(data, c('Hba-a1', 'Hba-a2', 'Hbb-bt', 'Hbb-bs'), ncol = 3)
dev.off()

# change back to working with peaks instead of gene activities
DefaultAssay(data) <- 'peaks'

Idents(data) <- 'predicted.id'
da_peaks <- FindMarkers(
  object = data,
  ident.1 = "Macrophages",
  ident.2 = "Dendritic cells",
  min.pct = 0.2,
  test.use = 'LR',
  latent.vars = 'peak_region_fragments'
)

head(da_peaks)

plot1 <- VlnPlot(
  object = data,
  features = rownames(da_peaks)[1],
  pt.size = 0.1,
  idents = c("Macrophages","Dendritic cells")
)
plot2 <- FeaturePlot(
  object = data,
  features = rownames(da_peaks)[1],
  pt.size = 0.1
)

plot1 | plot2

fc <- FoldChange(data, ident.1 = "Macrophages", ident.2 = "Dendritic cells")
head(fc)

open_mac <- rownames(da_peaks[da_peaks$avg_log2FC > 0.5, ])
open_dc <- rownames(da_peaks[da_peaks$avg_log2FC < -0.5, ])

closest_genes_mac <- ClosestFeature(data, regions = open_mac)
closest_genes_dc <- ClosestFeature(data, regions = open_dc)

head(closest_genes_mac)
head(closest_genes_dc)

CoveragePlot(
  object = data,
  region = rownames(da_peaks)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)

Idents(data) <- 'predicted.id'
pdf('UMAP_predicted_cell_types_atac.pdf', width = 6, height = 6)
DimPlot(
  object = data,
  label = F,
  repel = TRUE)
dev.off()

data_plot_1 <- data.frame(table(rna_reference$cell_type), 'Type' = 'RNA')
data_plot_1$Proportion <- data_plot_1$Freq/sum(data_plot_1$Freq)*100
data_plot_2 <- data.frame(table(data$predicted.id), 'Type' = 'ATAC')
data_plot_2$Proportion <- data_plot_2$Freq/sum(data_plot_2$Freq)*100
data_plot <- rbind(data_plot_1, data_plot_2)
colnames(data_plot)[1] <- 'Cell'
data_plot$Type <- factor(data_plot$Type, levels = c('RNA', 'ATAC'))

library(ggplot2)
library(ggpubr)
pdf('proportion.pdf', width = 5, height = 5)
ggbarplot(data_plot, x = "Cell", y = "Proportion",
          fill = "Type",
          width = 0.8,
          position = position_dodge(0.8),
          color = "white",
          palette = 'npg',
          label = round(data_plot$Proportion, 1),
          lab.vjust = 0.5,
          lab.hjust = 0.5
) + coord_flip() + theme(axis.text.y = element_text(hjust = 0)) + ylab('Proportion (%)') + xlab('')
dev.off()

save(data, file = 'annotated_atac.RData')

DefaultAssay(data) <- 'ACTIVITY'
pdf('dotplot_markers_final.pdf', width = 6, height = 12)
p <- DotPlot(data, features = c('Sftpa1', 'Sftpb', 'Sftpc', 'Sftpd',
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
p + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 0), axis.text.y = element_text(face = 'italic', hjust = 0), legend.position = 'bottom', legend.direction = 'vertical') + xlab('') + ylab('')
dev.off()

# convert to CellDataSet format and make the cicero object
data.cds <- as.cell_data_set(x = data)
data.cicero <- make_cicero_cds(data.cds, reduced_coordinates = reducedDims(data.cds)$UMAP)

# get the chromosome sizes from the Seurat object
genome <- seqlengths(data)

# convert chromosome sizes to a dataframe
genome.df <- data.frame("chr" = names(genome), "length" = genome)

# run cicero
conns <- run_cicero(data.cicero, genomic_coords = genome.df, sample_num = 100)
ccans <- generate_ccans(conns)
links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(data) <- links

save(conns, ccans, links, file = 'cicero.RData')
