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

##########
# Load RDS
##########
sce_QcCellsGenes_singlets = readRDS(file=snakemake@input[['rds_sce_cells_genes_singlets']]) # load data from RDS file
logcounts(sce_QcCellsGenes_singlets) = as.matrix(log2(counts(sce_QcCellsGenes_singlets) + 1)) # create logcounts assay from counts assay

#####################################
# Convert SCE object to Seurat object
#####################################
seurat_obj = as.Seurat(sce_QcCellsGenes_singlets, counts="counts", data="logcounts") # convert SCE to Seurat

#######################
# Lognorm normalization
#######################
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)

####################
# Save Seurat object
####################
saveRDS(seurat_obj,snakemake@output[['seurat_cells_genes_singlets_normed']])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])

#sce_QcCellsGenes_singlets = readRDS(file='OUTPUT/objects/sce/sce_cells_genes_singlets.rds')
#saveRDS(seurat_obj,"OUTPUT/objects/seurat/seurat_cells_genes_singlets_seurat_lognorm.rds")
#file.create("OUTPUT/.completion/step6_seurat_lognorm")