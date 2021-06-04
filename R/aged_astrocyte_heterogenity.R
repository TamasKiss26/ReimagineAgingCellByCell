library(Seurat)
library(tidyverse)
library(MAST)


## read data
BI_asc <- base::readRDS(file = './data/BI.ASC.normalized.reduced.seurat.rda')
BI_asc_aged <- base::subset(x = BI_asc, subset = age == '21-22mo')

BI_asc_aged <- Seurat::FindNeighbors(object = BI_asc_aged, dims = 1:20)
BI_asc_aged <- FindClusters(object = BI_asc_aged, resolution = 0.4)

BI_asc_aged_dp <- Seurat::DimPlot(object = BI_asc_aged) +
  ggplot2::labs(
    subtitle = 'Isolated from 21 to 22 Month Old Mouse Brain',
    title = 'Astrocytes'
  ) +
  ggplot2::xlab('UMAP 1') +
  ggplot2::ylab('UMAP 2') +
  ggplot2::theme(
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5)
  ) +
  ggplot2::coord_fixed(ratio = 1)

BI_asc_aged_dp

ggplot2::ggsave(
  filename = 'ASC_aged_dimplot.tiff',
  plot = BI_asc_aged_dp,
  device = 'tiff',
  scale = .9,
  dpi = 'retina'
)

BI_asc_aged_markers <- Seurat::FindAllMarkers(object = BI_asc_aged, test.use = 'MAST', only.pos = T)

BI_asc_aged_markers_top <- BI_asc_aged_markers %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::top_n(-2, p_val) %>% 
  dplyr::pull(gene)

Seurat::FeaturePlot(
  object = BI_asc_aged,
  features = BI_asc_aged_markers_top[c(1:7,9,11)]
)
