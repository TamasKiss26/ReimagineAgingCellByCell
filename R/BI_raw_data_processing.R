library(Seurat)
library(tidyverse)

# get expression data
count.raw <- readr::read_delim(
  file =  './data/raw/broad_institute_aging_mouse_brain_single_cell/expression_Aging_mouse_brain_portal_data_updated.txt', 
  delim = '\t'
)
count.raw <- count.raw%>%
  tibble::column_to_rownames('GENE')


# get metadata
meta_df <- readr::read_tsv(
  file = './data/raw/broad_institute_aging_mouse_brain_single_cell/meta_Aging_mouse_brain_portal_data.txt'
)

meta_df <- meta_df%>%
  dplyr::filter(NAME != 'TYPE')%>%
  tibble::column_to_rownames('NAME')%>%
  dplyr::select(cell_type, all_cells_by_age, nGene, nUMI)

base::colnames(meta_df) <- c('seurat_clusters', 'age','nGENE', 'nUMI') 

meta_df$seurat_clusters <- base::factor(meta_df$seurat_clusters)
meta_df$age             <- base::factor(meta_df$age)
meta_df$nGENE           <- base::as.numeric(meta_df$nGENE)
meta_df$nUMI            <- base::as.numeric(meta_df$nUMI)


# assemble Seurat object
BI <- Seurat::CreateSeuratObject(
  counts = count.raw,
  meta.data = meta_df,
  project = 'BI',
  min.cells = 0,
  min.features = 0,
  names.delim = '_',
  names.field = 6
)
BI <- Seurat::SetIdent(BI, value = 'seurat_clusters')

base::rm(count.raw, meta_df)

#BI <- readRDS('./data/BI.all.raw.count.seurat.rda')


# subset cells by cell type
cell_types <- Idents(BI) %>% levels()

BI_by_cell_types <- purrr::map(
  .x = cell_types,
  .f = function(ct){
    res <- base::subset(BI, cells = Seurat::WhichCells(object = BI, idents = ct))
    return(res)
  }
)
base::names(BI_by_cell_types) <- cell_types


# save raw count data
base::saveRDS(
  object = BI,
  file = './data/BI.all.raw.count.seurat.rda'
)

purrr::map(
  .x = cell_types,
  .f = function(ct){
    base::saveRDS(
      object = BI_by_cell_types[[ct]],
      file = paste0('./data/BI.', ct, '.raw.count.seurat.rda')
    )
  }
)


# normalize and reduce data
BI <- Seurat::NormalizeData(BI, normalization.method = "LogNormalize", scale.factor = 10000)
BI <- Seurat::FindVariableFeatures(BI, selection.method = "vst", nfeatures = 2000)
BI <- Seurat::ScaleData(BI)
BI <- Seurat::RunPCA(BI, features = VariableFeatures(object = BI))
BI <- Seurat::RunUMAP(BI, dims = 1:20)

cell_types <- cell_types[!cell_types %in% c("HypEPC", "NEUT", "TNC")]
BI_by_cell_types <-  purrr::map(
  .x = cell_types,
  .f = function(ct){
    base::print(ct)
    ssc <- BI_by_cell_types[[ct]]
    ssc <- Seurat::NormalizeData(ssc, normalization.method = "LogNormalize", scale.factor = 10000)
    ssc <- Seurat::FindVariableFeatures(ssc, selection.method = "vst", nfeatures = 2000)
    ssc <- Seurat::ScaleData(ssc)
    ssc <- Seurat::RunPCA(ssc, features = VariableFeatures(object = ssc))
    ssc <- Seurat::RunUMAP(ssc, dims = 1:20)
    return(ssc)
  }
)

names(BI_by_cell_types) <- cell_types

# save normalized and reduced data
base::saveRDS(
  object = BI,
  file = './data/BI.all.normalized.reduced.seurat.rda'
)

purrr::map(
  .x = cell_types,
  .f = function(ct){
    base::saveRDS(
      object = BI_by_cell_types[[ct]],
      file = paste0('./data/BI.', ct, '.normalized.reduced.seurat.rda')
    )
  }
)
