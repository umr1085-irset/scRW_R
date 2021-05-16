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
	counts <- assays(sce_QcCellsGenes_singlets)$counts # extract counts from SCE object
	rownames(counts) <-  rowData(sce_QcCellsGenes_singlets)$Symbol # extract rownames (cell barcodes)
	seurat_obj <- CreateSeuratObject(counts = counts, project = "scrw") # create Seurat object
	for (colname in colnames(colData(sce_QcCellsGenes_singlets))){ # loop over metadata columns in SCE object
		seurat_obj <- AddMetaData(object=seurat_obj, metadata=sce_QcCellsGenes_singlets[[colname]], col.name=colname) # add column
	}
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
        s.genes <- cc.genes$s.genes # extract genes associated to S cycle
        g2m.genes <- cc.genes$g2m.genes # extract genes associated to G2M cycle
        seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, assay = 'RNA', set.ident=TRUE) # compute cell cyle scores for all cells
	}
	seurat_obj <- FindVariableFeatures(seurat_obj) #FindVariableFeatures.
}

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

pdf(snakemake@output[['plot_pca_samples']]) # create PDF plot
DimPlot(seurat_obj, reduction = "pca", group.by = "Sample") + ggtitle("Normalized data colored by sample") # plot UMAP with samples
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

pdf(snakemake@output[['plot_umap_samples']]) # create PDF plot
DimPlot(seurat_obj, reduction = "umap", group.by = "Sample") + ggtitle("Normalized data colored by sample") # plot UMAP with samples
dev.off()

datalog = cbind(seurat_obj[["umap"]]@cell.embeddings, seurat_obj@meta.data[c("Sample","Phase","clust")])
write.csv(datalog, snakemake@output[["datalog_df"]]) # save pca coords to df for future report

####################
# Save Seurat object
####################
saveRDS(seurat_obj,file=snakemake@output[['rds_seurat']])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])

#pdf('OUTPUT/normalization/seurat_sctransform/umap_samples_normed_seurat_sctransform.pdf')
#DimPlot(seurat_obj, reduction = "umap", group.by = "Sample") + ggtitle("Normalized data colored by sample") # plot UMAP with cell phase
#dev.off()

#pdf('OUTPUT/normalization/seurat_sctransform/pca_samples_normed_seurat_sctransform.pdf') # create PDF plot
#DimPlot(seurat_obj, reduction = "pca", group.by = "Sample") + ggtitle("Normalized data colored by sample") # plot UMAP with samples
#dev.off()

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################
#
#
##################
## Import packages
##################
#library(Seurat)
#library(SingleCellExperiment)
#library(jackstraw)
#library(ggplot2)
#
###########
## Load RDS
###########
#if (snakemake@params[['seuratinput']]==0) { # if not a seurat input
#	sce_QcCellsGenes_singlets = readRDS(file=snakemake@input[['rds_sce_cells_genes_singlets_normed']]) # load data from RDS file
#	counts <- assays(sce_QcCellsGenes_singlets)$counts # extract counts from SCE object
#	rownames(counts) <-  rowData(sce_QcCellsGenes_singlets)$Symbol # extract rownames (cell barcodes)
#	seurat_obj <- CreateSeuratObject(counts = counts, project = "scrw") # create Seurat object
#	for (colname in colnames(colData(sce_QcCellsGenes_singlets))){ # loop over metadata columns in SCE object
#		seurat_obj <- AddMetaData(object=seurat_obj, metadata=sce_QcCellsGenes_singlets[[colname]], col.name=colname) # add column
#	}
#} else{
#	seurat_obj <- readRDS(file='OUTPUT/objects/seurat/seurat_cells_genes_singlets_seurat_sctransform.rds') # load data from RDS file
#}
#
##############################
## Scale and variable features
##############################
#if (snakemake@params[['scaled']]==0) {
#	if(snakemake@params[['regressoncellcyles']]){
#		s.genes <- cc.genes$s.genes # extract genes associated to S cycle
#		g2m.genes <- cc.genes$g2m.genes # extract genes associated to G2M cycle
#		seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, assay = 'RNA', set.ident=TRUE) # compute cell cyle scores for all cells
#		seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c('subsets_Mt_percent', "S.Score", "G2M.Score")) #ScaleData with no regression on cell cyles
#	} else {
#		seurat_obj <- ScaleData(seurat_obj,  vars.to.regress = c('subsets_Mt_percent')) #ScaleData with no regression on cell cyles
#	}
#	seurat_obj <- FindVariableFeatures(seurat_obj) #FindVariableFeatures.
#}
#
#############################
## Cell cycle phase detection
#############################
#
###########################
## Dimensionality reduction
###########################
#seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs=100)
#
#if (snakemake@params[['use_jack_straw']]==1) { # if JackStraw to be used
#	seurat_obj <- JackStraw(seurat_obj, reduction = "pca", num.replicate = 100, dims=100)
#	seurat_obj <- ScoreJackStraw(seurat_obj, dims = 1:100)
#	nDims <- which(seurat_obj[["pca"]]@jackstraw$overall.p.values[,2] > 0.001 )[1]-1
#	seurat_obj@reductions$pca@cell.embeddings <- seurat_obj@reductions$pca@cell.embeddings[,1:nDims]
#} else { # else if Seurat object with SCTransform normalization
#	nDims = 30
#}
#
#seurat_obj <- RunUMAP(seurat_obj, dims = 1:nDims) # run UMAP on PCA
#seurat_obj <- FindNeighbors(seurat_obj, dims = 1:nDims) # build knn graph then snn graph
#seurat_obj <- FindClusters(seurat_obj) # Louvain clustering
#clust.data <- sprintf("c%s",as.numeric(Idents(seurat_obj)))
#names(clust.data) <- names(Idents(seurat_obj))
#Idents(seurat_obj) <- clust.data
#seurat_obj@meta.data$clust <- Idents(seurat_obj)[colnames(seurat_obj)]
#
############
## Plot UMAP
############
#pdf(snakemake@output[['plot_pca_clusters']]) # create PDF plot
#DimPlot(seurat_obj, reduction="pca", label = TRUE) + NoLegend() # plot PCA with clusters
#dev.off()
#
#pdf(snakemake@output[['plot_pca_cellphase']]) # create PDF plot
#DimPlot(seurat_obj, reduction = "pca", group.by = "Phase") + ggtitle("Normalized data colored by cell phase") # plot PCA with cell phase
#dev.off()
#
############
## Plot UMAP
############
#pdf(snakemake@output[['plot_umap_clusters']]) # create PDF plot
#DimPlot(seurat_obj, reduction = "umap", label = TRUE) + NoLegend() # plot UMAP with clusters
#dev.off()
#
#pdf(snakemake@output[['plot_umap_cellphase']]) # create PDF plot
#DimPlot(seurat_obj, reduction = "umap", group.by = "Phase") + ggtitle("Normalized data colored by cell phase") # plot UMAP with cell phase
#dev.off()
#
#datalog = cbind(seurat_obj[["umap"]]@cell.embeddings, seurat_obj@meta.data[c("Sample","Phase","clust")])
#write.csv(datalog, snakemake@output[["datalog_df"]]) # save pca coords to df for future report
#
#####################
## Save Seurat object
#####################
#saveRDS(seurat_obj,file=snakemake@output[['rds_seurat']])
#
################
## Complete step
################
#file.create(snakemake@output[["step_complete"]])