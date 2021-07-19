library(dplyr)
library(Seurat)

###############################################
#  read in 10X data and create Seurat object  #
###############################################

cell_path<-'aggr_samples/outs/filtered_feature_bc_matrix'
cell_aggr<-Read10X(data.dir = cell_path)
seurat.object <- CreateSeuratObject(counts = cell_aggr, min.cells = 5)


#########################################
#   add metadata provided in aggr.csv   #
#########################################

# read in aggr.csv file used in cellranger aggr
aggr_csv<-read.csv("aggr.csv")

# get all cell barcodes
# these shoukld have a suffix like -1 or -2
# suffix of cell ID corresponds to which sample it came from in aggr.csv used to run cellranger
all_cells<-WhichCells(seurat.object)
# start an empty darta frame to deposit metadata into
metadata <- data.frame()

# loop through each row/suffix number
for (i in 1:nrow(aggr_csv)){
  suffix = paste("-", i, sep="")
  # count how many cell have that suffix number
  num_cells<-length(all_cells[endsWith(all_cells, suffix)])
  # take the row corresponding to that line number/sample and repeat by the number of cells
  cell_lines<-aggr_csv[rep(i, each = num_cells), ]
  # add these lines to the metadata dataframe
  metadata <- rbind(metadata, cell_lines)
}

# to be added to metadata, rownames of the dataframe need to be named after the cell barcode
rownames(metadata)<-all_cells

# add metadata to seurat object
seurat.object<-AddMetaData(seurat.object, metadata = metadata, col.name = colnames(metadata))


#########################################
#          Quality filtering            #
#########################################

# need to fill in format of Mt RNA (MT-/mt-/Mt-)
# assign percent MT RNA, use this and nFeatureRNA to filter the cells
seurat.object[["percent.mt"]]<-PercentageFeatureSet(object = seurat.object, pattern = "^MT-")
seurat.object<-subset(seurat.object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)


#########################################
#     Clustering and visualisation      #
#########################################

# Normalize, scale, cluster and visualise the data
seurat.object <- NormalizeData(seurat.object, verbose = FALSE)
seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000)
seurat.object <- ScaleData(seurat.object, verbose = FALSE)
seurat.object <- RunPCA(seurat.object, npcs = 10, verbose = FALSE)
seurat.object <- FindNeighbors(seurat.object, dims = 1:10)

# The resolution parameter in FindClusters sets the cluster granularity
# variable name changed to seurat.object_temp so that can change the resolution without having to rerun above code
seurat.object_temp <- FindClusters(seurat.object, resolution = 0.7)
seurat.object_temp <- RunTSNE(seurat.object_temp, dims = 1:10, check_duplicates = FALSE)


#########################################
#               Plotting                #
#########################################

# plot the seurat clusters
DimPlot(seurat.object_temp, reduction = "tsne", pt.size=0.7, order = TRUE)
