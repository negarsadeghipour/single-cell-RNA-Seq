# script to integrate scRNA-Seq datasets to correct for batch effects
# setwd("C:/Users/sade0010/OneDrive - University of Oklahoma/University of Oklahoma/bioinformatics practice/bioinformagician-scRNAseq/integration")

# load libraries
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# get data location

wd = getwd()
dirs <- list.dirs(path = wd, recursive = F, full.names = F)
dirs

# uncompressed the file

for(i in dirs){
  untar(paste0(i, '/', i))
  
}

dirs2 <- list.dirs(path = paste0(wd, '/data'), recursive = F, full.names = F)


for(x in dirs2){
  name <- gsub('_filtered_feature_bc_matrix', '', x)
  
  cnts <- ReadMtx(mtx = paste0(wd, '/data/', x, '/matrix.mtx.gz'),
          features = paste0(wd, '/data/', x, '/features.tsv.gz'),
          cells = paste0(wd, '/data/', x, '/barcodes.tsv.gz'))
  
  assign(name, CreateSeuratObject(counts = cnts))
}

# merge datasets

merged_seurat <- merge(HB17_background, y = c(HB17_PDX, HB17_tumor, HB30_PDX, HB30_tumor, 
                             HB53_background, HB53_tumor),
      add.cell.ids = ls()[4:10],
      project = 'HB')

merged_seurat



# QC and filtering --------------------------------------------------------

view(merged_seurat@meta.data)

# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)

merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Patient', 'Type', 'Barcode'),
         sep = '_')

# calculate mitochondrial percentage
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern = '^MT-')

# explore QC

# filtering
merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
         nFeature_RNA > 500 &
         mitoPercent < 10)

merged_seurat_filtered

# perform standard workflow steps to figure out if we see any batch effects
merged_seurat_filtered <- NormalizeData(object = merged_seurat_filtered)
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)

# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Patient')
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'Type',
              cols = c('red', 'green', 'blue'))
grid.arrange(p1, p2, ncol = 2, nrow = 2)

# perform integrating to correct for batch effects
obj.list <- SplitObject(merged_seurat_filtered,split.by = 'Patient')
for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}       

# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors
anchors <- FindIntegrationAnchors(object.list = obj.list,
                       anchor.features = features)

# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)

# scale dta, run PCA and UMAP and visualizze integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)

seurat.integrated <- RunPCA(object = seurat.integrated)

seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)

p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Patient')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'Type')

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
