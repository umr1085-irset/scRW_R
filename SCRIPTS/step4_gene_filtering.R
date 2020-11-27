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
#library(SingleCellExperiment)
library(scater)
library(scran)

################
# Load RDS files
################
sce = readRDS(file=snakemake@input[['rds_sce']])
sce_QcCells = readRDS(file=snakemake@input[['rds_sce_cells']])

###################################
# Filtering out low abundance genes
###################################
### Alternative approach: select genes that have non-zero counts in at least n cellules
print('Filtering out low abundance genes')
numcells <- nexprs(sce_QcCells, byrow=TRUE)
alt.keep <- numcells >= 10
sce_QcCellsGenes <- sce_QcCells[alt.keep,]

##########################
# Gene expression QC plots
##########################
print('Generating gene QC plots')
ave.counts <- rowMeans(counts(sce_QcCells)) # compute average count for each feature
pdf(snakemake@output[['QC_genes_AveCounts_PCAOutliersRemoved']]) # create PDF plot
hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
abline(v=log10(1), col="blue", lwd=2, lty=2)
dev.off()

# Number of expressing cells plot
pdf(snakemake@output[['QC_genes_NumCells_Counts_PCAOutliersRemoved']]) # create PDF plot
smoothScatter(log10(ave.counts), numcells, xlab=expression(Log[10]~"average count"),
              ylab="Number of expressing cells")
dev.off()

#Before cell PCA outlier removal:
ave_b <- calculateAverage(sce)
rowData(sce)$AveCount <- ave_b
pdf(snakemake@output[['QC_genes_AveGeneExpr_breaks100']]) # create PDF plot
hist(log10(ave_b), breaks=100, col="grey80")
dev.off()

#After cell outlier removal BUT before gene removal:
ave <- calculateAverage(sce_QcCells)
rowData(sce_QcCells)$AveCount <- ave
pdf(snakemake@output[['QC_genes_AveGeneExpr_PCAOutliersRemoved_breaks100']]) # create PDF plot
hist(log10(ave), breaks = 100, col="grey80")
dev.off()

#After cell outlier removal AND after gene removal:
ave <- calculateAverage(sce_QcCellsGenes)
rowData(sce_QcCellsGenes)$AveCount <- ave
pdf(snakemake@output[['QC_genes_AveGeneExpr_PCAOutliersRemoved_LowAbundanceGenesRemoved_breaks100']]) # create PDF plot
hist(log10(ave), breaks = 100, col="grey80")
dev.off()

# Detect the 50 most highly expressed genes:
# After cell and gene oulier removal
pdf(snakemake@output[['QC_genes_Top50Expr_GreyHist_QcCellsGenes']]) # create PDF plot
par(mar=c(5,4,1,1))
od1 = order(rowData(sce_QcCellsGenes)$AveCount, decreasing = TRUE)
barplot(rowData(sce_QcCellsGenes)$AveCount[od1[50:1]], las=1, 
        names.arg=rowData(sce_QcCellsGenes)$GeneName[od1[50:1]], 
        horiz=TRUE, cex.names=0.5, cex.axis=0.7, 
        xlab="ave # of UMI")
dev.off()

##########
# Save RDS
##########
saveRDS(sce_QcCellsGenes,snakemake@output[["rds_sce_cells_genes"]])
write.csv(rowData(sce_QcCellsGenes)[c('AveCount')], snakemake@output[["datalog_df"]])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])