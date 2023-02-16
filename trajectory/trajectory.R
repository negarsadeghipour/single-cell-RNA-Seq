# script to perform trajectory analysis 
# setwd("C:/Users/sade0010/OneDrive - University of Oklahoma/University of Oklahoma/bioinformatics practice/bioinformagician-scRNAseq/trajectory")
# data from: http://scrna.sklehabc.com/
# www.nature.com/articles/s41467-019-10291-0

set.seed(12324)

library(monocle3)
library(SeuratWrappers)
library(Seurat)
library(ggplot2)
library(tidyverse)

# read in data
markers <- read.delim('ABC_Marker.txt', header = T) # gene metadata

metadata <- read.delim('ABC_Meta.txt', header = T) # cell metadata

expr <- read.delim('ABC_umi_matrix_7551_cells.csv', header = T, sep = ',') # expression matrix


# creat a Seurat object -----------------------------------------------------------

expr.t <- t(expr)

seu.obj <- CreateSeuratObject(counts = expr.t)
View(seu.obj@meta.data)
# adding cell meta data information into this meta data
seu.obj@meta.data <- merge(seu.obj@meta.data, metadata, by.x = 'row.names', by.y = 'cell_id')
view(seu.obj@meta.data)
seu.obj@meta.data <- seu.obj@meta.data %>%
  column_to_rownames(var = 'Row.names')
seu.obj$mitopercent <- PercentageFeatureSet(seu.obj, pattern = '^MT-')  
seu.obj.filtered <- subset(seu.obj, subset = nCount_RNA > 800 &
                    nFeature_RNA > 500 &
                    mitopercent < 10)


# subset my seurat object to only contail B cells

unique(seu.obj.filtered@meta.data$population)

# changing seurat identity to population
Idents(seu.obj.filtered) <- seu.obj.filtered$population
b.seu <- subset(seu.obj.filtered, idents = "b")
b.seu
unique(b.seu@meta.data$redefined_cluster)

# pre-processing the data using seurat
b.seu <- NormalizeData(object = b.seu)
b.seu <- FindVariableFeatures(object = b.seu)
b.seu <- ScaleData(object = b.seu)
b.seu <- RunPCA(object = b.seu)
ElbowPlot(b.seu)
b.seu <- FindNeighbors(object = b.seu, dims = 1:30)
b.seu <- FindClusters(object = b.seu, resolution = 0.9)
b.seu <- RunUMAP(object = b.seu, dims = 1:30, n.neighbors = 50)

a1 <- DimPlot(b.seu, reduction = 'umap', group.by = 'redefined_cluster', label = T)
a2 <- DimPlot(b.seu, reduction = 'umap', group.by = 'seurat_clusters', label = T)
a1|a2


# mococle3 ----------------------------------------------------------------
# mococle3 requires cell_data_Set object
# convert seurat object to cell_data_set object for monocle3

cds <- as.cell_data_set(b.seu)
cds

# to get cell metadata
colData(cds)
# to gene metadata
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))

# to get counts
counts(cds)


# cluseter cells using clustering info from seurat's UMAP -----------------

# assign partition
recreate.partition <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)


cds@clusters$UMAP$partitions <- recreate.partition

# assign UMAP cpprdinate - cell embeddings
cds@int_colData@listData$reducedDims$UMAP <- b.seu@reductions$umap@cell.embeddings


# plot
cluster.before.trajectory <- plot_cells(cds,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
           color_cells_by = 'redefined_cluster',
           label_groups_by_cluster = FALSE,
           group_label_size = 5) +
  scale_color_manual(values = c('red', 'green', 'blue', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")

cluster.before.trajectory|cluster.names


# Learn trajectory graph --------------------------------------------------

cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'redefined_cluster',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5)


# order cells -------------------------------------------------------------

order_cells(cds, reduction_method = 'UMAP', root_cells = )
