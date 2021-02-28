library(Seurat)
library(tidyverse)
library(viridis)
library(tidymodels)

##===========##
## all cells ##
##===========##

# read data
BI <- base::readRDS(file = './data/BI.all.raw.count.seurat.rda')

# QC plots
Seurat::VlnPlot(object = BI, features = 'nCount_RNA') + Seurat::NoLegend()
Seurat::VlnPlot(object = BI, features = 'nFeature_RNA') + Seurat::NoLegend()
Seurat::VlnPlot(object = BI, features = 'nGENE') + Seurat::NoLegend()
Seurat::VlnPlot(object = BI, features = 'nUMI') + Seurat::NoLegend()

# scale and run dimension reduction on all cells
BI <- Seurat::FindVariableFeatures(BI, selection.method = "vst", nfeatures = 2000)
BI <- Seurat::ScaleData(BI)
BI <- Seurat::RunPCA(BI, features = VariableFeatures(object = BI))
BI <- Seurat::RunUMAP(BI, dims = 1:20)

# identify and plot microglia cells
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


fp_itgam <- Seurat::FeaturePlot(object = BI, features = 'Itgam', pt.size = .00001)+
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  ggtitle('Ox42 (Itgam) Expression') +
  scale_colour_viridis(direction = -1, option = 'viridis') +
  theme(
    legend.position = c(.9,.5)
  ) +
  coord_fixed(ratio = 1)

ggsave(
  filename = 'FeaturePlot_itgam.tiff',
  plot = fp_itgam,
  device = 'tiff',
  width = 7,
  height = 7,
  dpi = 'retina'
)


# identify and plot galectin
vp_lgals1 <- VlnPlot(object = BI, features = 'Lgals1', split.by = 'age') +
  ggtitle('Galectin-1 (Lgals1) Expression') +
  theme(
    axis.title.x = element_blank()
  )

ggsave(
  filename = 'VlnPlot_lgals1.tiff',
  plot = vp_lgals1,
  device = 'tiff',
  width = 9,
  height = 4,
  dpi = 'retina'
)


fp_lgals1 <- Seurat::FeaturePlot(object = BI, features = 'Lgals1', pt.size = .00001)+
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  ggtitle('Lgals1 (Galectin 1) Expression') +
  scale_colour_viridis(direction = -1, option = 'viridis') +
  theme(
    legend.position = c(.9,.5)
  ) +
  coord_fixed(ratio = 1)

ggsave(
  filename = 'FeaturePlot_lgals1.tiff',
  plot = fp_lgals1,
  device = 'tiff',
  width = 7,
  height = 7,
  dpi = 'retina'
)


##==========================##
## sub-populations of cells ##
##==========================##

# read data
mg.BI <- base::readRDS(file = './data/BI.MG.raw.count.seurat.rda' )

# calculate QC metrics
mg.BI[["percent.mt"]] <- PercentageFeatureSet(object = mg.BI, pattern = "^mt-")
mg.BI[["percent.rp"]] <- PercentageFeatureSet(object = mg.BI, pattern = "^Rp")

# scale and run PCA
mg.BI <- Seurat::FindVariableFeatures(object = mg.BI, selection.method = "vst", nfeatures = 2000)
mg.BI <- Seurat::ScaleData(object = mg.BI)
mg.BI <- Seurat::RunPCA(object = mg.BI, features = VariableFeatures(object = mg.BI))

# determine dimension
mg.BI <- Seurat::JackStraw(object = mg.BI, num.replicate = 100)
mg.BI <- Seurat::ScoreJackStraw(object = mg.BI, dims = 1:20)
Seurat::JackStrawPlot(object = mg.BI, dims = 1:20)

# UMAP and TSNE
mg.BI <- Seurat::RunTSNE(object = mg.BI, dims = 1:15)
mg.BI <- Seurat::RunUMAP(object = mg.BI, dims = 1:15)

# clustering
mg.BI <- Seurat::FindNeighbors(object = mg.BI, dims = 1:15)
mg.BI <- Seurat::FindClusters(object = mg.BI, resolution = 0.3)

mg.BI <- subset(mg.BI, cells = Seurat::WhichCells(object = mg.BI, idents = c(0:6)))

# plot QC metrics by cluster
Seurat::VlnPlot(object = mg.BI, features = 'percent.mt')
Seurat::VlnPlot(object = mg.BI, features = 'percent.rp')
Seurat::VlnPlot(object = mg.BI, features = 'nCount_RNA')
Seurat::VlnPlot(object = mg.BI, features = 'nFeature_RNA')

# DE genes
de_genes_mg <- Seurat::FindAllMarkers(object = mg.BI)

write.csv(file = 'de_genes_mg.csv', x = de_genes_mg)

# plot microglia subset
dp.mg <- DimPlot(object = mg.BI, reduction = 'umap', split.by = 'age') +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  ggtitle('Microglia Clusters') +
  theme(
    plot.title = element_text(hjust = .5),
    legend.position = c(.95,.6)
  ) +
  coord_fixed(ratio = 1)

dp.mg

ggsave(
  filename = 'DimPlot_microglia.tiff',
  plot = dp.mg,
  device = 'tiff',
  width = 7,
  height = 4,
  dpi = 'retina'
)

# count cells in the different populations
microglia_cluster_cell_n <- mg.BI[[]] %>%
  dplyr::select(seurat_clusters, age) %>%
  dplyr::group_by(seurat_clusters, age) %>%
  dplyr::summarise(n = n()) %>% 
  tidyr::pivot_wider(names_from = age, values_from = n)

write.csv(file = 'microglia_cluster_cell_n.csv', x = microglia_cluster_cell_n)

microglia_cluster_cell_n_plot <- microglia_cluster_cell_n %>%
  pivot_longer(!seurat_clusters, names_to = "age", values_to = "N") %>%
  ggplot(aes(x = seurat_clusters, y = N)) +
  geom_bar(aes(fill = age), position = "dodge", stat = "identity") +
  xlab('Microglia Clusters') +
  ylab('Number of Cells') +
  ggtitle('Microglia Clusters') +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.line = element_line(),
    legend.position = c(.8,.7),
    legend.text = element_text(size = 13),
    legend.title = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = .5, vjust = 1, face = 'bold', size = 15)
  )

microglia_cluster_cell_n_plot

ggsave(
  filename = 'microglia_cluster_cell_n.tiff',
  plot = microglia_cluster_cell_n_plot,
  device = 'tiff',
  width = 6,
  height = 4,
  dpi = 'retina'
)

# plot galectin positive cells
dp.mg_lgals1 <- DimPlot(
  object = mg.BI,
  cells.highlight = WhichCells(object = mg.BI, expression = Lgals1 > 0),
  split.by = 'age'
) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  ggtitle('Microglia Galectin 1 Expression')+
  scale_color_manual(
    label = c('Galectin 1 -', 'Galectin 1 +'),
    values = c('dark gray', 'red')
  ) +
  theme(
    plot.title = element_text(hjust = .5),
    legend.position = c(.88,.3)
  ) +
  coord_fixed(ratio = 1)

dp.mg_lgals1 

ggsave(
  filename = 'DimPlot_microglia_lgals1.tiff',
  plot = dp.mg_lgals1,
  device = 'tiff',
  width = 7,
  height = 4,
  dpi = 'retina'
)

# count galectin positive cells in the different populations
galectin_poz <- WhichCells(object = mg.BI, expression = Lgals1 > 0)
mg.BI[['galectin']] <- ifelse(colnames(mg.BI) %in% galectin_poz, 'poz', 'neg')

microglia_cluster_cell_galectin_positive_n_plot <- mg.BI[[]] %>%
  dplyr::select(seurat_clusters, galectin) %>%
  dplyr::group_by(seurat_clusters, galectin) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = galectin, values_from = n) %>%
  dplyr::transmute(
    poz_ratio = poz / (neg + poz) * 100,
    neg_ratio = neg / (neg + poz) * 100
    ) %>%
  tidyr::pivot_longer(!seurat_clusters, names_to = 'galectin', values_to = 'ratio') %>%
  ggplot(aes(x = seurat_clusters, y = ratio, fill = galectin)) +
  geom_bar(position = 'fill', stat = "identity") +
  xlab('Microglia Clusters') +
  ylab('Ratio of Gelectin 1 Positive Cells') +
  ggtitle('Microglia Galectin 1 Expression by Cluster') +
  scale_fill_manual(
    label = c('Galectin 1 -', 'Galectin 1 +'),
    values = c('dark gray', 'red')
  ) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.line = element_line(),
    legend.position = c(.8,.7),
    legend.text = element_text(size = 13),
    legend.title = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(hjust = .5, vjust = 1, face = 'bold', size = 15)
  )

microglia_cluster_cell_galectin_positive_n_plot 

ggsave(
  filename = 'microglia_cluster_gaelctin_pos.tiff',
  plot = microglia_cluster_cell_galectin_positive_n_plot ,
  device = 'tiff',
  width = 6,
  height = 4,
  dpi = 'retina'
)


# model galectin positivity
galectin_train <- mg.BI[[]] %>% 
  dplyr::select(seurat_clusters, galectin) %>%
  dplyr::mutate_if(
    .predicate = is.character,
    .funs = factor
  )

glm_fit <- logistic_reg() %>%
  set_engine('glm') %>%
  fit(galectin ~ seurat_clusters, data = galectin_train)

tidy(glm_fit)

new_clusters <- tibble(
  seurat_clusters = unique(galectin_train$seurat_clusters)
)

mean_pred <- predict(
  glm_fit,
  new_data = new_clusters,
  type = 'prob'
)

conf_int <- predict(
  glm_fit,
  new_data = new_clusters,
  type = 'conf_int'
)

cluster_result <- new_clusters %>%
  bind_cols(mean_pred) %>%
  bind_cols(conf_int)







#=============================================================================================
fp_lgals3 <- FeaturePlot(object = mg.BI, features = 'Lgals3', reduction = 'umap') + 
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  ggtitle('Lgals3 Expression in Microglia') +
  scale_colour_viridis(direction = -1, option = 'viridis') +
  theme(
    plot.title = element_text(hjust = .5)
  ) +
  coord_fixed(ratio = 1)


fp_tnf <- FeaturePlot(object = mg.BI, features = 'Tnf', reduction = 'umap') + 
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  ggtitle('Tnf Expression in Microglia') +
  scale_colour_viridis(direction = -1, option = 'viridis') +
  theme(
    plot.title = element_text(hjust = .5)
  ) +
  coord_fixed(ratio = 1)

ggsave(
  filename = 'FeaturePlot_microglia_tnf.tiff',
  plot = fp_tnf,
  device = 'tiff',
  width = 8,
  height = 4,
  dpi = 'retina'
)



tnf_poz <- WhichCells(object = mg.BI, expression = Tnf > 0)
mg.BI[['Tnf']] <- ifelse(colnames(mg.BI) %in% tnf_poz, 'poz', 'neg')

mg.BI[[]] %>%
  dplyr::select(seurat_clusters, Tnf) %>%
  dplyr::group_by(seurat_clusters, Tnf) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = Tnf, values_from = n) %>%
  dplyr::mutate(ratio = poz/(neg+poz)*100) %>%
  ggplot(aes(x = seurat_clusters, y = ratio)) +
  geom_bar(stat = "identity") +
  xlab('Microglia Cluster') +
  ylab('Ratio of Tnf Positive Cells') +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.line = element_line(),
    legend.position = c(.9,.7),
    legend.text = element_text(size = 13),
    legend.title = element_blank(),
    panel.background = element_blank()
  )


Seurat::DimPlot(
  object = mg.BI,
  cells.highlight = Seurat::WhichCells(object = mg.BI, expression = Nrp1 > 0)
)

Seurat::FeaturePlot(
  object = mg.BI,
  features = 'Nrp1'
)

Nrp1_poz <- WhichCells(object = mg.BI, expression = Nrp1 > 0)
mg.BI[['Nrp1']] <- ifelse(colnames(mg.BI) %in% Nrp1_poz, 'poz', 'neg')

mg.BI[[]] %>%
  dplyr::select(seurat_clusters, Nrp1) %>%
  dplyr::group_by(seurat_clusters, Nrp1) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = Nrp1, values_from = n) %>%
  dplyr::mutate(ratio = poz/(neg+poz)*100) %>%
  ggplot(aes(x = seurat_clusters, y = ratio)) +
  geom_bar(stat = "identity") +
  xlab('Microglia Cluster') +
  ylab('Ratio of Nrp1 Positive Cells') +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.line = element_line(),
    legend.position = c(.9,.7),
    legend.text = element_text(size = 13),
    legend.title = element_blank(),
    panel.background = element_blank()
  )










FeaturePlot(object = mg.BI, features = 'Lgals1', reduction = 'umap') +
  scale_color_viridis(direction = -1)
FeaturePlot(object = mg.BI, features = 'Tnf', reduction = 'umap') +
  scale_color_viridis(direction = -1)

FeaturePlot(object = mg.BI, features = 'Tnf', reduction = 'umap') +
  scale_color_viridis(direction = -1)

FeaturePlot(object = mg.BI, features = 'Lgals1', reduction = 'umap')
FeaturePlot(object = mg.BI, features = 'Lgals1', split.by = 'age')


library(GENIE3)
library(RcisTarget)
library(AUCell)

mg.BI.matrix <- Seurat::GetAssayData(mg.BI)
mg.BI.matrix <- as.matrix(mg.BI.matrix)
weightMat <- GENIE3(mg.BI.matrix)

##

VlnPlot(object = mg.BI, features = 'Lgals1', split.by = 'age') +
  ggtitle('Galectin-1 (Lgals1) Expression') +
  theme(
    axis.title.x = element_blank()
  )

Seurat::FeaturePlot(object = mg.BI, features = 'Tnf', pt.size = .00001)+
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  ggtitle('Tnf Expression') +
  scale_colour_viridis(direction = -1, option = 'viridis') +
  theme(
    legend.position = c(.9,.5)
  ) +
  coord_fixed(ratio = 1)


##





##


pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)



galectin_poz <- WhichCells(object = BI, expression = Lgals1 > 0)
BI[['galectin_poz']] <- ifelse(colnames(BI) %in% galectin_poz, 'poz', 'neg')

galectin_poz_tab <- BI[[]] %>%
  dplyr::select(galectin_poz, age, seurat_clusters) %>%
  dplyr::group_by(age, seurat_clusters) %>%
  dplyr::summarise(n = n())



mg <- readRDS('./data/aging.MG.seurat.rda')

DimPlot(object = mg)
FeaturePlot(object = mg, features = 'Lgals1', split = 'age')


mg[[]] %>%
  colnames()

galectin_poz <- WhichCells(object = mg, expression = Lgals1 > 0)
mg[['galectin']] <- ifelse(colnames(mg) %in% galectin_poz, 'poz', 'neg')

FeaturePlot(object = mg, features = 'Lgals1')
DimPlot(object = mg)
mg.BI[[]] %>%
  dplyr::select(seurat_clusters, galectin) %>%
  dplyr::group_by(seurat_clusters, galectin) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = galectin, values_from = n) %>%
  dplyr::mutate(ratio = poz/(neg+poz)*100)

DimPlot(object = mg, split = 'age')
mg[[]] %>%
  dplyr::select(seurat_clusters, age) %>%
  dplyr::group_by(seurat_clusters, age) %>%
  dplyr::summarise(n = n()) %>%
  tidyr::pivot_wider(names_from = age, values_from = n)


DEgenes <- Seurat::FindAllMarkers(object = mg)
DEgenes %>% 
  dplyr::filter(cluster == 2) %>%
  View()
