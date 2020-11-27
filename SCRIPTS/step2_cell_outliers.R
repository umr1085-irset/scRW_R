#############################################################
# PREPROCESSING SCRIPT
# P. RIVAUD
# 2020/04
#############################################################

# In R, Snakefile parameters can be accessed with the following:
# print(snakemake@input[[1]])
# print(snakemake@output[["rds_sce"]])
# file.create(snakemake@output[["rds_sce"]])

###############
# Load packages
###############
library(scater)

###############
# Load RDS file
###############
sce = readRDS(file=snakemake@input[['rds_sce']])

#######################
# PCA outlier detection
#######################
print('Removing outlier cells')
sce  <- runColDataPCA(sce,variables=c("subsets_Ri_percent","subsets_Mt_percent"),outliers=TRUE,name='PCA_coldata') # run pca on col data
pcacoords <- as.data.frame(reducedDim(sce, "PCA_coldata")[,1:2])
pcacoords$outlier <- sce$outlier

pdf(snakemake@output[['outliers_plot']]) # create PDF plot
plotReducedDim(sce, dimred="PCA_coldata", colour_by="outlier") # plot PCA
dev.off()

sce_QcCells <- sce[,which(sce$outlier==FALSE)] # keep non-outlier cells
# saveRDS(sce,snakemake@input[['rds_sce']]) # overwrite sce object file
saveRDS(sce_QcCells,snakemake@output[['rds_sce_cells']]) # save cell-filtered sce object
write.csv(pcacoords, snakemake@output[["datalog_df"]]) # save pca coords to df for future report

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])

#write.csv(pcacoords,'OUTPUT/report/datalog/datalogstep2_df.csv')