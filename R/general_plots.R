library(Seurat)
library(tidyverse)

# read data
BI <- readRDS('./data/BI.all.normalized.reduced.seurat.rda')


# combined DimPlot
dp <- Seurat::DimPlot(object = BI, label = T, label.box = T, label.size = 3, repel = T, shuffle = T) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  ggtitle('Cell Type Clusters') +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = .5)
  ) +
  coord_fixed(ratio = 1)

ggsave(
  filename = 'DimPlot.tiff',
  plot = dp,
  device = 'tiff',
  width = 7,
  height = 7,
  dpi = 'retina'
)
