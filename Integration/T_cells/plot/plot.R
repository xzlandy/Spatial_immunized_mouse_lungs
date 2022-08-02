library(Seurat)
library(Signac)

setwd('~/Box/RWorkSpace/Spatial_combine/Integration/T_cells/plot/')

load('../tcell.RData')

DefaultAssay(data) <- 'peaks'

# load('../../cicero.RData')
# Links(data) <- links

pdf('peaks_Il6.pdf', width = 5, height = 5)
CoveragePlot(
  object = data,
  region = 'Il6',
  extend.upstream = 2000,
  extend.downstream = 2000
)
dev.off()

pdf('peaks_Il6ra.pdf', width = 5, height = 5)
CoveragePlot(
  object = data,
  region = 'Il6ra',
  extend.upstream = 2000,
  extend.downstream = 2000
)
dev.off()

pdf('peaks_Cxcl5.pdf', width = 5, height = 5)
CoveragePlot(
  object = data,
  region = 'Cxcl5',
  extend.upstream = 4000,
  extend.downstream = 4000
)
dev.off()

pdf('peaks_Tgfbr1.pdf', width = 5, height = 5)
CoveragePlot(
  object = data,
  region = 'Tgfbr1',
  extend.upstream = 2000,
  extend.downstream = 2000
)
dev.off()

pdf('peaks_Tgfbr2.pdf', width = 5, height = 5)
CoveragePlot(
  object = data,
  region = 'Tgfbr2',
  extend.upstream = 2000,
  extend.downstream = 2000
)
dev.off()

pdf('peaks_Acvr1b.pdf', width = 5, height = 5)
CoveragePlot(
  object = data,
  region = 'Acvr1b',
  extend.upstream = 2000,
  extend.downstream = 2000
)
dev.off()