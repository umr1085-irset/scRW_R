#############################################################
# PREPROCESSING SCRIPT
# P. RIVAUD
# 2020/04
#############################################################

# In R, Snakefile parameters can be accessed with the following:
# print(snakemake@input[[1]])
# print(snakemake@output[["rds_sce"]])
# file.create(snakemake@output[["rds_sce"]])

#################
# Import packages
#################
library(Seurat)
library(SingleCellExperiment)
library(jackstraw)

##########
# Load RDS
##########
sce_QcCellsGenes_singlets = readRDS(file=snakemake@input[['rds_sce_cells_genes_singlets']]) # load data from RDS file
logcounts(sce_QcCellsGenes_singlets) = as.matrix(log2(counts(sce_QcCellsGenes_singlets) + 1)) # create logcounts assay from counts assay
seurat_obj <- as.Seurat(sce_QcCellsGenes_singlets, counts="counts", data="logcounts")

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs=100)
seurat_obj <- JackStraw(seurat_obj, reduction = "pca", num.replicate = 100, dims=100)
seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:100)
nDims <- which( seurat_obj[["pca"]]@jackstraw$overall.p.values[,2] > 0.001 )[1]-1

seurat_obj@reductions$pca@cell.embeddings <- seurat_obj@reductions$pca@cell.embeddings[,1:nDims]
seurat_obj <- RunUMAP(seurat_obj, dims = 1:nDims)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:nDims)
seurat_obj <- FindClusters(seurat_obj)

clust.data <- sprintf("c%s",as.numeric(Idents(seurat_obj)))
names(clust.data) <- names(Idents(seurat_obj))

Idents(seurat_obj) <- clust.data
seurat_obj@meta.data$clust <- Idents(seurat_obj)[colnames(seurat_obj)]

####################
# Save Seurat object
####################
saveRDS(seurat_obj,file=snakemake@output[['rds_seurat']])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])

#marker.data <- FindAllMarkers(subset.seurat_obj, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

# R TEST -------------------------------------------------------------------------------------------------
#path_ = 'OUTPUT/objects/sce/sce_cells_genes_singlets.rds'
#sce_QcCellsGenes_singlets = readRDS(file=path_)
#sce_QcCellsGenes_singlets = runPCA(sce_QcCellsGenes_singlets) # library(scater)
#logcounts(sce_QcCellsGenes_singlets) = as(log2(counts(sce_QcCellsGenes_singlets) + 1), "sparseMatrix") # create logcounts assay from counts assay
#logcounts(sce_QcCellsGenes_singlets) = as.matrix(log2(counts(sce_QcCellsGenes_singlets) + 1)) # create logcounts assay from counts assay
#seurat_obj <- as.Seurat(sce_QcCellsGenes_singlets, counts="counts", data="logcounts")