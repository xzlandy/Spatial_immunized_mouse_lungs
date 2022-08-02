library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(openxlsx)
library(stringr)
library(DOSE)

setwd('~/Box/RWorkSpace/Spatial_combine/DE_analysis/')

# SET THE DESIRED ORGANISM HERE
library(org.Mm.eg.db)

# reading in data from deseq2
gsea_plot <- function(path){
  df <- read.xlsx(path)
  df <- df[df$avg_log2FC != 0,]
  write.xlsx(df, path, rowNames = T)
  
  # we want the log2 fold change 
  original_gene_list <- df$avg_log2FC
  
  # name the vector
  names(original_gene_list) <- df$gene
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  gse <- gseGO(geneList=gene_list, 
               ont ="BP", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = org.Mm.eg.db, 
               pAdjustMethod = "fdr")
  
  # require(DOSE)
  # dotplot(gse, showCategory=10, split=".sign", color = "pvalue") + facet_grid(.~.sign)
  
  gse_m <- gse@result
  gse_m <- gse_m[order(gse_m$NES),]
  gse_m_filt = rbind(head(gse_m, n = 30),
                     tail(gse_m, n = 30 ))
  gse_m_filt <- gse_m_filt[order(gse_m_filt$NES, decreasing = F),]
  gse_m_filt$Description <- factor(gse_m_filt$Description, levels = gse_m_filt$Description)
  g = ggplot(gse_m_filt, aes(Description, NES)) +
    geom_segment(aes(xend=Description, y=0, yend=NES)) +
    geom_point(aes( fill = p.adjust, size = setSize),
               shape=21, stroke=2) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title="GSEA - Biological Processes") + 
    theme_minimal()
  return(g)
}


files <- list.files(pattern = 'mast')

for(i in files){
  pdf(paste0('gseGO_', str_split_fixed(i, '_DE_mast.xlsx', 2)[1,1], '.pdf'), width = 10, height = 10)
  print(gsea_plot(i))
  dev.off()
}

# emapplot(gse, showCategory = 10)
# 
# # categorySize can be either 'pvalue' or 'geneNum'
# cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)
# 
# ridgeplot(gse) + labs(x = "enrichment distribution")
# 
# # Use the `Gene Set` param for the index in the title, and as the value for geneSetId
# gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)