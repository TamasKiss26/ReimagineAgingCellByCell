library(Seurat)
library(tidyverse)
library(viridis)

BI <- readRDS('./data/BI.all.normalized.reduced.seurat.rda')


FeaturePlot(
  object = BI,
  features = 'Wdr62'
)

id <- str_detect(
  string = rownames(BI),
  pattern = 'Wdr6'
) 

rownames(BI)[id]



FeaturePlot(
  object = BI,
  features = 'Cdk5rap2'
) + 
  scale_colour_viridis(direction = -1, option = 'viridis')

