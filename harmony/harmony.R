# script to integrate across conditions using Harmony 
# setwd("C:/Users/sade0010/OneDrive - University of Oklahoma/University of Oklahoma/bioinformatics practice/bioinformagician-scRNAseq/harmony")

# set seed for reproducibility
set.seed(1234)

library(harmony)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(SeuratData)


# get data ----------------------------------------------------------------

# to see which datasets are available
AvailableData()

# install dataset
InstallData("ifnb")

# load dataset
LoadData("ifnb")

str(ifnb)

# QC and filtering

ifnb$mito.percent <- PercentageFeatureSet(ifnb, pattern = '^MT-')
view(ifnb@meta.data)
# explore QC
# filter
ifnb.filtered <- subset(ifnb, subset = nCount_RNA > 800 &
                          nFeature_RNA > 200 &
                          mito.percent < 5)

# standard workflow steps
ifnb.filtered <- NormalizeData(ifnb.filtered)
ifnb.filtered <- FindVariableFeatures(ifnb.filtered)
ifnb.filtered <- ScaleData(ifnb.filtered)
ifnb.filtered <- RunPCA(ifnb.filtered)
ElbowPlot(ifnb.filtered)
ifnb.filtered <- RunUMAP(ifnb.filtered, dims = 1:20, reduction = 'pca')

before <- DimPlot(ifnb.filtered, reduction = 'umap', group.by = 'stim')

# run Harmony
ifnb.harmony <- ifnb.filtered %>%
  RunHarmony(group.by.vars = 'stim', plot_convergance = FALSE)

ifnb.harmony@reductions

ifnb.harmony.embed <- Embeddings(ifnb.harmony, "harmony")
ifnb.harmony.embed[1:10, 1:10]

# UMAP and clustering
ifnb.harmony <- ifnb.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# visualize
after <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim')

before|after

view(ifnb.harmony@meta.data)

clusters <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
condition <- DimPlot(ifnb.harmony, reduction = 'umap', group.by = 'stim')

condition|clusters

# findAll markers
# better for data where one group is present -- this is the represenation of how it works on this data
FindAllMarkers(ifnb_harmony,
               logfc.threshold = 0.25,
               min.pct = 0.1,
               only.pos = TRUE,
               test.use = 'DESeq2',
               slot = 'counts')

# findConserved markers

DefaultAssay(ifnb.harmony)
markers_cluster3 <- FindConservedMarkers(ifnb.harmony,
                     ident.1 = 4,
                     grouping.var = 'stim')

head(markers_cluster3)

# visualize top features
FeaturePlot(ifnb.harmony, features = c('FCGR3A'), min.cutoff = 'q10')

# min-cut off explanation
seq(1,5)
SetQuantile('q50', seq(1,5))
SetQuantile('q10', seq(1,5))

# rename quantile 3 ident
Idents(ifnb.harmony)
ifnb.harmony <- RenameIdents(ifnb.harmony, `4` = 'CD16 Mono')

DimPlot(ifnb.harmony, reduction = 'umap', label = TRUE)

# cells already annotated in the meta data !!! not a real world senario
Idents(ifnb.harmony) <- ifnb.harmony@meta.data$seurat_annotations
Idents(ifnb.harmony)

DimPlot(ifnb.harmony, reduction = 'umap', label = TRUE)

# findMarkers between conditions
ifnb.harmony$celltype.cnd <- paste0(ifnb.harmony$seurat_annotations, '_', ifnb.harmony$stim)
view(ifnb.harmony)
Idents(ifnb.harmony) <- ifnb.harmony$celltype.cnd
DimPlot(ifnb.harmony, reduction = 'umap', label = T)

# find markers
n.interferon.response <- FindMarkers(ifnb.harmony, ident.1 = 'CD16 Mono_STIM', ident.2 = 'CD16 Mono_CTRL')

head(n.interferon.response)
# plotting conserved features vs DE features between conditions
head(markers_cluster3)

FeaturePlot(ifnb.harmony, features = c('FCGR3A', 'AIF1', 'IFIT1'), split.by = 'stim', min.cutoff = 'q10')
