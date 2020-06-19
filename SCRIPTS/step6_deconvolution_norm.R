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
library(scran)
library(scater)

##########
# Load RDS
##########
sce_QcCellsGenes_filtered = readRDS(file=snakemake@input[['rds_sce_cells_genes_filtered']]) # load data from RDS file

##################################
# Subset data from sample singlets
##################################
if (snakemake@params[['individual_samples']]) { # if individual samples mode, subset
	singlets = read.table(snakemake@input[['sample_singlets']]) # read list of barcodes
	singlets = as.character(singlets$V1) # to vector
	sce_QcCellsGenes_filtered = sce_QcCellsGenes_filtered[,singlets] # subset cells
}

####################
# Perform clustering
####################
clusters = quickCluster(sce_QcCellsGenes_filtered, min.mean=0.1, method="igraph") # quickCluster: cluster cells based on rank correlations in their gene expression profiles
saveRDS(clusters,snakemake@output[["rds_clusters"]])

#####################
# Compute sum factors
#####################
sce_QcCellsGenes_filtered = computeSumFactors(sce_QcCellsGenes_filtered, cluster=clusters, min.mean=0.1)
sce_QcCellsGenes_filtered_normed = logNormCounts(sce_QcCellsGenes_filtered)

pdf(snakemake@output[['plot_Norm_HistSizeFactors']])
hist(sizeFactors(sce_QcCellsGenes_filtered), breaks=100, col="grey80")
dev.off()

pdf(snakemake@output[['plot_Norm_SizeFactorsVsTotalCountsPerMillion']])
plot(sce_QcCellsGenes_filtered$total/1e6, sizeFactors(sce_QcCellsGenes_filtered),log="xy", xlab="Library size (millions)", ylab="Size factor", cex=0.3, pch=20, col=rgb(0.1,0.2,0.7,0.3))
dev.off()

pdf(snakemake@output[['plot_Norm_SizeFactorsVsTotalCounts_smooth']]) 
smoothScatter(sce_QcCellsGenes_filtered$total, sizeFactors(sce_QcCellsGenes_filtered), xlab="total counts", ylab="size factors", cex=0.3, pch=20, col=rgb(0.1,0.2,0.7,0.3))
dev.off()

#############
# Save to RDS
#############
saveRDS(sce_QcCellsGenes_filtered_normed, snakemake@output[['rds_sce_cells_genes_filtered_normed']])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])