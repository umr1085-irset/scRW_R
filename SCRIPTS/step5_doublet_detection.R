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
library(DoubletFinder)
library(SingleCellExperiment)
library(Seurat)

########################################
# Extract individual sample barcode list
########################################
bcsgz=gzfile(snakemake@input[['bcsfile']],'rt') # barcodes file handle
bcs=read.csv(bcsgz,header=F) # read barcodes
bcs=bcs$V1 # extract barcodes to vector
idx=match(snakemake@wildcards[["sample"]],snakemake@params[['samplelist']]) # get index of sample in aggregation sample list
bcs=gsub("1", paste0(idx), bcs) # replace default "1" index in individual barcode list with sample index from aggregation list

##########
# Load RDS
##########
sce_QcCellsGenes = readRDS(file=snakemake@input[['rds_sce_cells_genes']]) # load data from RDS file
logcounts(sce_QcCellsGenes) = as.matrix(log2(counts(sce_QcCellsGenes) + 1)) # create logcounts assay from counts assay
#logcounts(sce_QcCellsGenes) = log2(counts(sce_QcCellsGenes) + 1) # create logcounts assay to convert to Seurat object

###########################################
# Extract cells matching the current sample
###########################################
allbcs=rownames(colData(sce_QcCellsGenes)) # get all barcodes from sce_QcCellsGenes object
intersection=intersect(allbcs,bcs) # get intersection between original sample barcode list and all barcodes from sce_QcCellsGenes (after cell filtering)
sub = sce_QcCellsGenes[,intersection] # subset sce_QcCellsGenes object using cell barcodes from current sample
sub.seurat <- as.Seurat(sub, counts="counts", data="logcounts")

#################################
# Preparation steps DoubletFinder
#################################
sub.seurat <- NormalizeData(sub.seurat)
sub.seurat <- ScaleData(sub.seurat)
sub.seurat <- FindVariableFeatures(sub.seurat, selection.method = "vst", nfeatures = 2000)
sub.seurat <- RunPCA(sub.seurat)
#sub.seurat <- FindClusters(sub.seurat, dims.use = 1:30, resolution = 0.5, print.output = FALSE)
#conda install -c conda-forge umap-learn
#Error in RunUMAP.default(object = data.use, assay = assay, n.neighbors = n.neighbors,  : 
#  Cannot find UMAP, please install through pip (e.g. pip install umap-learn).
#sub.seurat <- RunUMAP(sub.seurat, dims = 1:10)
sweep.res.list_sub.seurat <- paramSweep_v3(sub.seurat, PCs = 1:10, sct = FALSE)
sweep.stats_sub.seurat<- summarizeSweep(sweep.res.list_sub.seurat, GT = FALSE)
#BCmvn = mean-variance-normalized bimodality coefficient
bcmvn_sub.seurat <- find.pK(sweep.stats_sub.seurat)
#DoubletFinder is sensitive to heterotypic doublets -- i.e., doublets formed from transcriptionally-distinct cell states -- 
# but is insensitive to homotypic doublets -- i.e., doublets formed from transcriptionally-similar cell states.
#BCmvn distributions also feature a single maximum for scRNA-seq datasets generated without sample-multiplexing, enabling pK selection.
#pk <- as.numeric(as.vector(bcmvn_sub.seurat$pK)[which.max(bcmvn_sub.seurat$BCmetric)]) #0.005
pk = 0.09
nExp_poi <- round(0.075 * dim(sub.seurat)[2]) #191
#pN ~ This defines the number of generated artificial doublets, expressed as a proportion of the merged real-artificial data.
# Default is set to 25%, based on observation that DoubletFinder performance is largely pN-invariant

###################
# Run DoubletFinder
###################
sub.seurat <- doubletFinder_v3(sub.seurat, PCs = 1:10, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) # run DoubletFinder
indexDFInfo <- which(names(sub.seurat@meta.data) == sprintf("DF.classifications_0.25_%s_%s",pk,nExp_poi)) # extract DF results column index
DoubletsIndices <- which(sub.seurat@meta.data[indexDFInfo][,1] =="Doublet") # extract cell indices classified as doublets
SingletsIndices <- which(sub.seurat@meta.data[indexDFInfo][,1] =="Singlet") # extract cell indices classified as singlets
Doublets <- rownames(sub.seurat@meta.data[DoubletsIndices,]) # extract cell barcodes classified as doublets
Singlets <- rownames(sub.seurat@meta.data[SingletsIndices,]) # extract cell barcodes classified as singlets
write.table(Doublets, file=snakemake@output[['doublets_sample_file']], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) # save doublets barcodes to file
write.table(Singlets, file=snakemake@output[['singlets_sample_file']], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) # save singlets barcodes to file

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])