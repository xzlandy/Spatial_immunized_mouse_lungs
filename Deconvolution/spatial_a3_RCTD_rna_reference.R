library(RCTD)
library(Matrix)
library(Seurat)

setwd('~/Spatial/Deconvolution/')

puck <- readRDS('spatial_a3/STRef.rds')
barcodes <- colnames(puck@counts) #pixels to be used (a list of barcode names). 

reference = readRDS('rna_reference/SCRef.rds')

myRCTD <- create.RCTD(puck, reference, max_cores = 24, test_mode = F, CELL_MIN_INSTANCE = 6)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

results <- myRCTD@results
save(myRCTD, file = 'spatial_a3_RCTD_rna_reference.RData')
# normalize the cell type proportions to sum to 1.
norm_weights = as.matrix(sweep(results$weights, 1, rowSums(results$weights), '/'))
boxplot(norm_weights)

plot_frac = lapply(colnames(norm_weights), function(i) plot_puck_continuous(puck, barcodes, norm_weights[,i], title = i))

