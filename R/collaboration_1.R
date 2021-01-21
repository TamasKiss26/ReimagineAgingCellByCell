library(Seurat)
library(tidyverse)
library(viridis)
library(ggpubr)

BI <- readRDS('./data/BI.all.normalized.reduced.seurat.rda')

fp1 <- Seurat::FeaturePlot(
  object = BI,
  features = c('Cbs', 'Cth', 'Mpst'),
  split.by = 'age',
  pt.size = .1
)

ggplot2::ggsave(
  filename = 'FeaturePlot_Cbs_Cth_Mpst.tiff',
  plot = fp1,
  device = 'tiff',
  width = 12,
  height = 18,
  dpi = 'retina'
)

FindMarkers(
  object = BI,
  subset.ident = "ASC",
  ident.1 = '21-22mo',
  group.by = 'age'
)
