# script to perform standard workflow steps to analyze single cell RNA-Seq data
# data: 10k. Human PBMC Multiome v1.0, Chromium Controller

setwd("C:/Users/sade0010/OneDrive - University of Oklahoma/University of Oklahoma/bioinformatics practice/bioinformagician-scRNAseq/ssgeneexpression")


# load libraries ----------------------------------------------------------

library(Seurat)
library(tidyverse)


# read the counts data and load it into seurat object ---------------------

fname = dir()

nsclc.sparse.m <- Read10X_h5(filename = fname)

str(nsclc.sparse.m)

cts <- nsclc.sparse.m$`Gene Expression`



# Initialize the Seurat obect with the raw (non-normalized data) ----------

nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)

str(nsclc.seurat.obj)
nsclc.seurat.obj


# 1. QC ----------------------------------------------------------------------

view(nsclc.seurat.obj@meta.data)
# % MT reads
# add the percentage of mitochondrial genes as a column to the seurat object

nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
view(nsclc.seurat.obj@meta.data)

# see the genes as violin plot
VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_smooth(method = 'lm')


# 2. Filtering ------------------------------------------------------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                             percent.mt < 5)


# 3. Normalize data -------------------------------------------------------

nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj)


# 4. identify highly variable features ---------------------------------------

nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# identify the 10 most variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)

# plot variable features with and without labels

plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


# 5. Scaling --------------------------------------------------------------

all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)              



# 6. perform PCA analysis -------------------------------------------------

nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(nsclc.seurat.obj, dims = 1:2
           , cells = 500, balance = TRUE)

# determine dimensionality of the data
ElbowPlot(nsclc.seurat.obj)


# 7. clustering -----------------------------------------------------------

nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# understanding resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))
view(nsclc.seurat.obj@meta.data)


DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.1", label = TRUE)

# setting identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)



# non-linear dimensionality reduction -------------------------------------

nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
DimPlot(nsclc.seurat.obj, reduction = "umap")
