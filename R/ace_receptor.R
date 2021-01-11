library(Seurat)
library(tidyverse)
library(viridis)
library(ggpubr)

BI <- readRDS('./data/BI.all.normalized.reduced.seurat.rda')

# combined FeaturePlot
MyFeaturePlot <- function(gene_abr, gene_name){
  mfp <- Seurat::FeaturePlot(object = BI, features = gene_abr, pt.size = .00001) +
    xlab('UMAP 1') +
    ylab('UMAP 2') +
    ggtitle(gene_name) +
    scale_colour_viridis(direction = -1, option = 'viridis') +
    theme(
      legend.position = c(.9,.2)
    ) +
    coord_fixed(ratio = 1)
  
  ggsave(
    filename = paste0('FeaturePlot_', gene_abr, '.tiff'),
    plot = mfp,
    device = 'tiff',
    width = 7,
    height = 7,
    dpi = 'retina'
  )
}

MyFeaturePlot(gene_abr = 'Agtrap', gene_name = 'angiotensin II, type I receptor-associated protein')
MyFeaturePlot(gene_abr = 'Agt', gene_name = 'angiotensinogen')
MyFeaturePlot(gene_abr = 'Agtr1b', gene_name = 'angiotensin II receptor, type 1b')
MyFeaturePlot(gene_abr = 'Arap1', gene_name = 'ArfGAP with RhoGAP domain, ankyrin repeat and PH domain 1 ')
MyFeaturePlot(gene_abr = 'Arrb2', gene_name = 'arrestin, beta 2')


# split FeaturePlot
MySplitFeaturePlot <- function(gene_abr, gene_name){
  mfp_young <- Seurat::FeaturePlot(object = Seurat::SplitObject(object = BI, split.by = 'age')[[1]], features = gene_abr, pt.size = .00001) +
    xlab('UMAP 1') +
    ylab('UMAP 2') +
    ggtitle('2-3mo') +
    scale_colour_viridis(direction = -1, option = 'viridis') +
    theme(
      legend.position = c(.9,.2)
    ) +
    coord_fixed(ratio = 1)
  
  mfp_aged <- Seurat::FeaturePlot(object = Seurat::SplitObject(object = BI, split.by = 'age')[[2]], features = gene_abr, pt.size = .00001) +
    xlab('UMAP 1') +
    ylab('UMAP 2') +
    ggtitle('21-22mo') +
    scale_colour_viridis(direction = -1, option = 'viridis') +
    theme(
      legend.position = c(.9,.2)
    ) +
    coord_fixed(ratio = 1)
  
  mfp_comb <- ggpubr::ggarrange(
    plotlist = list(mfp_young, mfp_aged),
    ncol = 2, nrow = 1
  ) %>% annotate_figure(p = ., top = text_grob(gene_name, face = 'bold', size = 18))
  
  ggsave(
    filename = paste0('SplitFeaturePlot_', gene_abr, '.tiff'),
    plot = mfp_comb,
    device = 'tiff',
    width = 14,
    height = 7.5,
    dpi = 'retina'
  )
}

MySplitFeaturePlot(gene_abr = 'Agtrap', gene_name = 'angiotensin II, type I receptor-associated protein')
MySplitFeaturePlot(gene_abr = 'Agt', gene_name = 'angiotensinogen')
MySplitFeaturePlot(gene_abr = 'Agtr1b', gene_name = 'angiotensin II receptor, type 1b')
MySplitFeaturePlot(gene_abr = 'Arap1', gene_name = 'ArfGAP with RhoGAP domain, ankyrin repeat and PH domain 1 ')
MySplitFeaturePlot(gene_abr = 'Arrb2', gene_name = 'arrestin, beta 2')


# dotplot
DP_young <- DotPlot(object = Seurat::SplitObject(object = BI, split.by = 'age')[[1]], features = c('Agtrap', 'Agt', 'Agtr1b', 'Arap1', 'Arrb2')) +
  xlab('') +
  ylab('') +
  NoLegend()

DP_aged <- DotPlot(object = Seurat::SplitObject(object = BI, split.by = 'age')[[2]], features = c('Agtrap', 'Agt', 'Agtr1b', 'Arap1', 'Arrb2')) +
  xlab('') +
  ylab('') +
  NoLegend()

ggpubr::get_legend(DP_aged)

ggpubr::ggarrange(
  plotlist = list(DP_young, DP_aged),
  ncol = 2, nrow = 1
)
