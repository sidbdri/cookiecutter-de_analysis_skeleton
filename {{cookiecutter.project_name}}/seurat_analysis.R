library(dplyr)
library(Seurat)

cell_path<-'aggr_samples/outs/filtered_feature_bc_matrix'
cell_aggr<-Read10X(data.dir = cell_path)

seurat.object <- CreateSeuratObject(counts = cell_aggr, min.cells = 5) # include features in at least 5 cells

# need to fill in format of Mt RNA (MT/mt/Mt)
seurat.object[["percent.mt"]]<-PercentageFeatureSet(object = seurat.object, pattern = "^MT-")

seurat.object<-subset(seurat.object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)
#seurat.object<-subset(seurat.object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
# Normnalise, scale, cluster and represent cells using tsne dim reduction
seurat.object <- NormalizeData(seurat.object, verbose = FALSE)
seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000)
seurat.object <- ScaleData(seurat.object, verbose = FALSE)
seurat.object <- RunPCA(seurat.object, npcs = 10, verbose = FALSE)
seurat.object <- FindNeighbors(seurat.object, dims = 1:10)


seurat.object_temp <- FindClusters(seurat.object, resolution = 0.7)
seurat.object_temp <- RunTSNE(seurat.object_temp, dims = 1:10, check_duplicates = FALSE)

DimPlot(seurat.object_temp, reduction = "tsne", pt.size=0.7, order = TRUE)
