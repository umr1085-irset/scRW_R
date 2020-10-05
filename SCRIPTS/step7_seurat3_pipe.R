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
if (snakemake@params[['seuratinput']]==0) { # if not a seurat input
	sce_QcCellsGenes_singlets = readRDS(file=snakemake@input[['rds_sce_cells_genes_singlets_normed']]) # load data from RDS file
	logcounts(sce_QcCellsGenes_singlets) = as.matrix(log2(counts(sce_QcCellsGenes_singlets) + 1)) # create logcounts assay from counts assay
	seurat_obj <- as.Seurat(sce_QcCellsGenes_singlets, counts="counts", data="logcounts")	
} else{
	seurat_obj <- readRDS(file=snakemake@input[['rds_sce_cells_genes_singlets_normed']]) # load data from RDS file
}

##########################
# Dimensionality reduction
##########################
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs=100) # run PCA

if (snakemake@params[['seuratinput']]==0) { # if not a seurat input
	seurat_obj <- JackStraw(seurat_obj, reduction = "pca", num.replicate = 100, dims=100)
	seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:100)
	nDims <- which(seurat_obj[["pca"]]@jackstraw$overall.p.values[,2] > 0.001 )[1]-1
	seurat_obj@reductions$pca@cell.embeddings <- seurat_obj@reductions$pca@cell.embeddings[,1:nDims]	
} else { # else if Seurat object with SCTransform normalization
	nDims = 30
}

seurat_obj <- RunUMAP(seurat_obj, dims = 1:nDims) # run UMAP on PCA obj
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:nDims) # 
seurat_obj <- FindClusters(seurat_obj)
clust.data <- sprintf("c%s",as.numeric(Idents(seurat_obj)))
names(clust.data) <- names(Idents(seurat_obj))
Idents(seurat_obj) <- clust.data
seurat_obj@meta.data$clust <- Idents(seurat_obj)[colnames(seurat_obj)]

###########
# Plot UMAP
###########
pdf(snakemake@output[['plot_umap']]) # create PDF plot
DimPlot(seurat_obj, label = TRUE) + NoLegend() # plot UMAP with clusters
dev.off()

####################
# Save Seurat object
####################
saveRDS(seurat_obj,file=snakemake@output[['rds_seurat']])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])