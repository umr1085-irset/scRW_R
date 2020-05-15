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

###################
# Run DoubletFinder
###################

print(snakemake@input[["sample"]])
ith_QcCellsGenes <- doubletFinder(ith_QcCellsGenes, PCs = 1:30, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE)

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])


####################################
### Doublet detection            ###
####################################
#https://github.com/chris-mcginnis-ucsf/DoubletFinder
# DoubletFinder can be broken up into 4 steps:
# (1) Generate artificial doublets from existing scRNA-seq data
# (2) Pre-process merged real-artificial data
# (3) Perform PCA and use the PC distance matrix to find each cell's proportion of artificial k nearest neighbors (pANN)
# (4) Rank order and threshold pANN values according to the expected number of doublets

###3- Ensure that input data is cleared of low-quality cell clusters AND then of low QC genes:
### 3.1- Reduce the sce object of the single sample to the cells quality controlled on the whole expression matrix (7 samples at once) (from script step4_ITH_ccRCC_20190529_QCs.R)
sce_QcCells <- sce[,QcCells]
sprintf("Number of cells in this samples reduced from: %s to %s after quality control", length(colData(sce)$X),length(QcCells))
sprintf("%s cells were removed during the quality control on cells and genes of this sample", length(setdiff(colData(sce)$X,QcCells))) #290

### 3.2- Reduce the sce object to the genes quality controlled on the whole expression matrix (7 samples at once) ###
QC_genes <- rownames(sce_QcCellsGenes_ALL) #GeneNames
sce_QcCellsGenes <- sce_QcCells[QC_genes,]
###/3

### 4- Convert the sce object {scater} into a {Seurat} object ---
### Need to have the logcounts slot to create the Seurat object:
### Create an additional slot with log-transformed counts: =log2(DeepthNormalizedData!)
#NEED TO NAME THIS "logcounts" TO BE ABLE TO SIMPLY ACCESS THE DATA WITH logcounts(Selected.Genes.Cells.sce)
assay(sce_QcCellsGenes,"logcounts") <- log2(counts(sce_QcCellsGenes) + 1) # logcounts slot created.

### Creation of the ith (IntraTumoralHeterogeneity) Seurat object:
ith_QcCellsGenes <- Convert(from = sce_QcCellsGenes, to = "seurat")

### 5- DoubletFinder steps: inspired from: cf https://github.com/chris-mcginnis-ucsf/DoubletFinder
### 5.1- Pre-process Seurat object with standard method: ---
ith_QcCellsGenes <- NormalizeData(ith_QcCellsGenes)
ith_QcCellsGenes <- ScaleData(ith_QcCellsGenes)

### The redundant genes were removed in the QC step on cells ans genes
# RedundantGenes <- QC_genes[which(duplicated(QC_genes)=="TRUE")]
# RedundantGenesIndices <-which((rownames(ith_QcCellsGenes@data) %in% RedundantGenes)==TRUE)
# sprintf("Nb of unredundant genes: %s", nrow(ith_QcCellsGenes@data[-RedundantGenesIndices,]))
# ith_QcCellsGenes@data <- ith_QcCellsGenes@data[-RedundantGenesIndices,]
ith_QcCellsGenes <- FindVariableGenes(ith_QcCellsGenes, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = FALSE) #output saved in ith_QcCellsGenes@hvg.info

ith_QcCellsGenes <- RunPCA(ith_QcCellsGenes, pc.genes = ith_QcCellsGenes@var.genes, pcs.print = 0, pcs.compute = 30)
ith_QcCellsGenes <- FindClusters(ith_QcCellsGenes, dims.use = 1:30, resolution = 0.5, print.output = FALSE)
#conda install -c conda-forge umap-learn, ne fonctionne pas sur fenetre conda en parallele
ith_QcCellsGenes <- RunUMAP(ith_QcCellsGenes, dims=1:30) #warnings.warn(errors.NumbaDeprecationWarning(msg, self.func_ir.loc))

### 5.2- pK Identification (no ground-truth): ---
####remove gene redundancy in raw.data: No more because done at the previous QC on cells and genes step!
#ith_QcCellsGenes@raw.data <- ith_QcCellsGenes@raw.data[-RedundantGenesIndices,]
#inspired from: https://www.gitmemory.com/chris-mcginnis-ucsf
sweep.res.list_ith_QcCellsGenes <- paramSweep(ith_QcCellsGenes, PCs = 1:30)
sweep.stats_ith_QcCellsGenes    <- summarizeSweep(sweep.res.list_ith_QcCellsGenes, GT = FALSE)
#BCmvn = mean-variance-normalized bimodality coefficient
bcmvn_ith_QcCellsGenes <- find.pK(sweep.stats_ith_QcCellsGenes)

#DoubletFinder is sensitive to heterotypic doublets -- i.e., doublets formed from transcriptionally-distinct cell states -- 
# but is insensitive to homotypic doublets -- i.e., doublets formed from transcriptionally-similar cell states.
#BCmvn distributions also feature a single maximum for scRNA-seq datasets generated without sample-multiplexing, enabling pK selection.
pk <- as.numeric(as.vector(bcmvn_ith_QcCellsGenes$pK)[which.max(bcmvn_ith_QcCellsGenes$BCmetric)]) #0.005
nExp_poi <- round(0.075 * length(ith_QcCellsGenes@cell.names)) #191
#pN ~ This defines the number of generated artificial doublets, expressed as a proportion of the merged real-artificial data.
# Default is set to 25%, based on observation that DoubletFinder performance is largely pN-invariant

### 5.3- Run DoubletFinder
#pk should be adjusted for each datset
ith_QcCellsGenes <- doubletFinder(ith_QcCellsGenes, PCs = 1:30, pN = 0.25, pK = pk, nExp = nExp_poi, reuse.pANN = FALSE)
#output info in: ith_QcCellsGenes@meta.data$DF.classifications_0.25_0.005_191