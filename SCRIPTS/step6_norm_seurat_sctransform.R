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

##########
# Load RDS
##########
sce_QcCellsGenes_singlets = readRDS(file=snakemake@input[['rds_sce_cells_genes_singlets']]) # load data from RDS file
logcounts(sce_QcCellsGenes_singlets) = as(log2(counts(sce_QcCellsGenes_singlets) + 1), "sparseMatrix") # create logcounts assay from counts assay

#####################################
# Convert SCE object to Seurat object
#####################################
seurat_obj <- as.Seurat(sce_QcCellsGenes_singlets, counts="counts", data="logcounts")












###############################################################
### R TEST
###############################################################
path_ = 'OUTPUT/objects/sce/sce_cells_genes_singlets.rds'
sce_QcCellsGenes_singlets = readRDS(file=path_)
logcounts(sce_QcCellsGenes_singlets) = as(log2(counts(sce_QcCellsGenes_singlets) + 1), "sparseMatrix") # create logcounts assay from counts assay
seurat_obj <- as.Seurat(sce_QcCellsGenes_singlets, counts="counts", data="logcounts")

# run sctransform
#seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "subsets_Mt_percent", verbose = FALSE)

seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
saveRDS(seurat_obj,'seurat_obj.rds')
seurat_obj = readRDS('seurat_obj.rds')
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)
DimPlot(seurat_obj, label = TRUE) + NoLegend()

# problem
# which version of the data to use?
# using mt percent variable from step1

#> warnings()
#Messages d'avis :
#1: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#2: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#3: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#4: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#5: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#6: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#7: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#8: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#9: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#10: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#11: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#12: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#13: In sqrt(1/i) : production de NaN
#14: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#15: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#16: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#17: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#18: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#19: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#20: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#21: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#22: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#23: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#24: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#25: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#26: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#27: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#28: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#29: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#30: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#31: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#32: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#33: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#34: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#35: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#36: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#37: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#38: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#39: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#40: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#41: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#42: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#43: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#44: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#45: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#46: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#47: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#48: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#49: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached
#50: In theta.ml(y = y, mu = fit$fitted) : iteration limit reached

