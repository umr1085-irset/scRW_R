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
options(future.globals.maxSize = snakemake@params[["future_globals_maxsize"]] * 1024^2)
options(future.globals.maxSize = 2000 * 1024^2)
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
saveRDS(seurat_obj,'tmp_seurat_obj.rds')
###########################
# SCTransform normalization
###########################

## normalize data with SCTransform()
if(snakemake@params[['regressoncellcyles']]){ # if set to normalize regressing on cell cycles
	s.genes <- cc.genes$s.genes # extract genes associated to S cycle
	g2m.genes <- cc.genes$g2m.genes # extract genes associated to G2M cycle
	seurat_obj <- SCTransform(seurat_obj, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('subsets_Mt_percent'))	# normalize before computing cell cyle scores
	seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes, assay = 'SCT', set.ident=TRUE) # compute cell cyle scores for all cells
	seurat_obj <- SCTransform(seurat_obj, assay = 'RNA', new.assay.name = 'SCT', vars.to.regress = c('subsets_Mt_percent', 'S.Score', 'G2M.Score')) # normalize again but this time including also the cell cycle scores
} else {
	seurat_obj = SCTransform(seurat_obj, vars.to.regress = "subsets_Mt_percent", verbose = FALSE) # run Seurat sctransform method
}

####################
# Save Seurat object
####################
saveRDS(seurat_obj,snakemake@output[['seurat_cells_genes_singlets_normed']])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])