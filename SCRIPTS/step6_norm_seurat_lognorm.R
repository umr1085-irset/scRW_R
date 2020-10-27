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

#####################################
# Convert SCE object to Seurat object
#####################################
counts <- assays(sce_QcCellsGenes_singlets)$counts # extract counts from SCE object
rownames(counts) <-  rowData(sce_QcCellsGenes_singlets)$Symbol # extract rownames (cell barcodes)
seurat_obj <- CreateSeuratObject(counts = counts, project = "scrw") # create Seurat object
for (colname in colnames(colData(sce_QcCellsGenes_singlets))){ # loop over metadata columns in SCE object
	seurat_obj <- AddMetaData(object=seurat_obj, metadata=sce_QcCellsGenes_singlets[[colname]], col.name=colname) # add column
}

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