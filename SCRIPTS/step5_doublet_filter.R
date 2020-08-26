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
library(SingleCellExperiment)

##########
# Load RDS
##########
sce_QcCellsGenes = readRDS(file=snakemake@input[['rds_sce_cells_genes']]) # load data from RDS file

#######################################
# Retrieve singlets barcodes from files
#######################################
singlets=vector() # conta
for (f in snakemake@input[['bcs_file_list']]){ # for each file in input list (one file per sample)
	l=read.table(f) # read list of barcodes
	l=as.character(l$V1) # to vector
	singlets=c(singlets,l) # concatenate to vector
}
write.table(singlets, file=snakemake@output[['all_singlets']], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE) # save list of singlet barcodes to file

################################################
# Save sce object with DoubletFinder information
################################################
sce_QcCellsGenes_DF = sce_QcCellsGenes
colData(sce_QcCellsGenes_DF)["DoubletFinder"] <- "Doublet"
colData(sce_QcCellsGenes_DF)["DoubletFinder"][singlets,] <- "Singlet"
saveRDS(sce_QcCellsGenes_DF,snakemake@output[['rds_sce_cells_genes_DF']])

####################################
# Save sce object with singlets only
####################################
rds_sce_cells_genes_singlets = sce_QcCellsGenes_DF[,singlets]
saveRDS(rds_sce_cells_genes_singlets, snakemake@output[['rds_sce_cells_genes_singlets']])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])