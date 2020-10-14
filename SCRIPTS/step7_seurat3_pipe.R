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
library(ggplot2)

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

#############################
# Scale and variable features
#############################
if (snakemake@params[['scaled']]==0) {
	if(snakemake@params[['regressoncellcyles']]){
		s.genes <- cc.genes$s.genes # extract genes associated to S cycle
		g2m.genes <- cc.genes$g2m.genes # extract genes associated to G2M cycle
		seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, assay = 'RNA', set.ident=TRUE) # compute cell cyle scores for all cells
		seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c('subsets_Mt_percent', "S.Score", "G2M.Score")) #ScaleData with no regression on cell cyles
	} else {
		seurat_obj <- ScaleData(seurat_obj,  vars.to.regress = c('subsets_Mt_percent')) #ScaleData with no regression on cell cyles
	}
	seurat_obj <- FindVariableFeatures(seurat_obj) #FindVariableFeatures.
}

############################
# Cell cycle phase detection
############################

##########################
# Dimensionality reduction
##########################
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs=100)

if (snakemake@params[['use_jack_straw']]==1) { # if JackStraw to be used
	seurat_obj <- JackStraw(seurat_obj, reduction = "pca", num.replicate = 100, dims=100)
	seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:100)
	nDims <- which(seurat_obj[["pca"]]@jackstraw$overall.p.values[,2] > 0.001 )[1]-1
	seurat_obj@reductions$pca@cell.embeddings <- seurat_obj@reductions$pca@cell.embeddings[,1:nDims]
} else { # else if Seurat object with SCTransform normalization
	nDims = 30
}

seurat_obj <- RunUMAP(seurat_obj, dims = 1:nDims) # run UMAP on PCA
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:nDims) # build knn graph then snn graph
seurat_obj <- FindClusters(seurat_obj) # Louvain clustering
clust.data <- sprintf("c%s",as.numeric(Idents(seurat_obj)))
names(clust.data) <- names(Idents(seurat_obj))
Idents(seurat_obj) <- clust.data
seurat_obj@meta.data$clust <- Idents(seurat_obj)[colnames(seurat_obj)]

###########
# Plot UMAP
###########
pdf(snakemake@output[['plot_pca_clusters']]) # create PDF plot
DimPlot(seurat_obj, reduction="pca", label = TRUE) + NoLegend() # plot PCA with clusters
dev.off()

pdf(snakemake@output[['plot_pca_cellphase']]) # create PDF plot
DimPlot(seurat_obj, reduction = "pca", group.by = "Phase") + ggtitle("Normalized data colored by cell phase") # plot PCA with cell phase
dev.off()

###########
# Plot UMAP
###########
pdf(snakemake@output[['plot_umap_clusters']]) # create PDF plot
DimPlot(seurat_obj, reduction = "umap", label = TRUE) + NoLegend() # plot UMAP with clusters
dev.off()

pdf(snakemake@output[['plot_umap_cellphase']]) # create PDF plot
DimPlot(seurat_obj, reduction = "umap", group.by = "Phase") + ggtitle("Normalized data colored by cell phase") # plot UMAP with cell phase
dev.off()

####################
# Save Seurat object
####################
saveRDS(seurat_obj,file=snakemake@output[['rds_seurat']])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])