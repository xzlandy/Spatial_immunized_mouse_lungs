library(Signac)
library(Seurat)
library(Matrix)
library(ggplot2)
library(patchwork)

setwd('~/Box/RWorkSpace/Spatial_combine/Integration/T_cells/SCENIC/')

load('../transferred_rna_reference.RData')

exprMat <- as.matrix(GetAssayData(rna_reference, slot = 'counts', assay = 'RNA'))
cellInfo <- data.frame(seuratCluster=rna_reference$cell_type)

### Initialize settings
library(SCENIC)
library(AUCell)
scenicOptions <- initializeScenic(org="mgi", dbDir="../cisTarget_databases", nCores=40, dbs = c("mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather", "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"))
# scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)

### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)
scenicOptions@settings$nCores <- 1
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
scenicOptions@settings$nCores <- 40

# Optional: Binarize activity
# aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_log)
# savedSelections <- shiny::runApp(aucellApp)
# newThresholds <- savedSelections$thresholds
# scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
# saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

# Export:
# saveRDS(cellInfo, file=getDatasetInfo(scenicOptions, "cellInfo")) # Temporary, to add to loom
export2loom(scenicOptions, exprMat)

# To save the current status, or any changes in settings, save the object again:
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

### Exploring output
# Check files in folder 'output'
# Browse the output .loom file @ http://scope.aertslab.org

# output/Step2_MotifEnrichment_preview.html in detail/subset:
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Stat1"]
viewMotifs(tableSubset)

# output/Step2_regulonTargetsInfo.tsv in detail:
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Stat1" & highConfAnnot==TRUE]
viewMotifs(tableSubset)

# Cell-type specific regulators (RSS):
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "seuratCluster"], )
rssPlot <- plotRSS(rss)
pdf('SCENIC_regulator.pdf', width = 8, height = 32)
# plotly::ggplotly(rssPlot$plot)
rssPlot$plot
dev.off()

rss_subset <- rss[,c(7,15,2,3)]
rssPlot_subset <- plotRSS(rss_subset)
pdf('SCENIC_regulator_subset.pdf', width = 4, height = 16)
# plotly::ggplotly(rssPlot$plot)
rssPlot_subset$plot + theme(text = element_text(size = 16))
dev.off()

AUC=getAUC(regulonAUC)
rna_reference[['SCENIC']] <- CreateAssayObject(AUC)

rna_tcell <- subset(rna_reference, cell_type == 'Th17' | cell_type == 'Th1' | cell_type == 'Other T cells 1' | cell_type == 'Other T cells 2' )
rna_tcell$cell_type <- factor(as.character(rna_tcell$cell_type), levels = c('Th17', 'Th1', 'Other T cells 1', 'Other T cells 2'))

rna_tcell <- SCTransform(rna_tcell, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
# These are now standard steps in the Seurat workflow for visualization and clustering
rna_tcell <- RunPCA(rna_tcell, verbose = FALSE)
ElbowPlot(rna_tcell, ndims = 50)

rna_tcell <- RunUMAP(rna_tcell, dims = 1:50, verbose = FALSE)

rna_tcell <- FindNeighbors(rna_tcell, dims = 1:50, verbose = FALSE)
rna_tcell <- FindClusters(rna_tcell, resolution = 1)

DimPlot(rna_tcell, label = TRUE) + NoLegend()
DefaultAssay(rna_tcell) <- 'SCENIC'

pdf('scatter_Rora-extended.pdf', width = 8, height = 8)
FeaturePlot(rna_tcell, 'Rora-extended (15g)')
dev.off()

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

pdf('SCENIC_regulator_heatmap.pdf', width = 8, height = 40)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
dev.off()

