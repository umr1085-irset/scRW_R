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
#library(DoubletFinder)
library(SingleCellExperiment)
#library(Seurat)

##########
# Load RDS
##########
sce_QcCellsGenes = readRDS(file=snakemake@input[['rds_sce_cells_genes']]) # load data from RDS file
singlets=vector()
for (f in snakemake@input[['bcs_file_list']]) {
    l=read.table(f)
    l=as.character(l$V1)
    singlets=c(singlets,l)
}
print(singlets)
print(typeof(singlets))
write.table(singlets, file=snakemake@output[['all_singlets']], sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
sce_QcCellsGenes_filtered = sce_QcCellsGenes[,singlets]

# Save RDS
saveRDS(sce_QcCellsGenes_filtered,snakemake@output[['rds_sce_cells_genes_filtered']])

#### 4.2- Doublets prediction on each biological sample
#CONVERSIONFILE="/groups/irset/archives/SingleCell/projects/HumanRCC/ITH_ccRCC/201904/GenewizName2SampleName_20190409.txt"
#ConversionInfo <- as.matrix(read.table(CONVERSIONFILE, header=F, sep="\t"))
#SampleNames    <- ConversionInfo[,2]
#SampleIDs      <- matrix(unlist(strsplit(ConversionInfo[,1],"_")),nrow=length(SampleNames),ncol=5,byrow=T)[,1]
#
#sceData <- sce_QcCellsGenes_ALL
#for (i in 1:nrow(ConversionInfo)) {
#	SampleName <- SampleNames[i]
#	SampleID   <- SampleIDs[i]
#	cat("DoubletDetection ON:: SampleName:", SampleName, "and SampleID:", SampleID, "\n")
#	DoubletDetection(i,SampleName,SampleID,sceData)
#}
#
#### 4.3- Save the DoubletFinder information in a single file ###
#OUTFile <- "ITH_ccRCC_DoubletFinder.txt"
#for (SampleName in SampleNames) {
#	DFFile <- sprintf("ITH_ccRCC_DoubletFinder_%s.txt",SampleName)
#	Data <- as.matrix(read.table(DFFile))[,1]
#	sprintf("Sample: %s --> Number of doublets: %s.", SampleName, length(Data))
#
#	write.table(Data,file = OUTFile,sep = "\t", append=TRUE, quote= FALSE, row.names=FALSE, col.names=FALSE)
#}
#################################################################
#
#
#### 4.4- Save the DoubletFinder information in the sce object on all the samples ###
#ALL_Doublets <- as.matrix(read.table(OUTFile))[,1]
#length(ALL_Doublets) #2095
#
#sce_QcCellsGenes_DF <- sce_QcCellsGenes_ALL 
#sce_QcCellsGenes_DF # 21144 genes x 27921 cells, only counts slot!!
##Add these information in the sce object
#colData(sce_QcCellsGenes_DF)["DoubletFinder"] <- "Singlet"
#colData(sce_QcCellsGenes_DF)["DoubletFinder"][ALL_Doublets,] <- "Doublet"
#table(colData(sce_QcCellsGenes_DF)["DoubletFinder"])
##Doublet Singlet
##   2095   25826
#
##Save RDS:
#saveRDS(sce_QcCellsGenes_DF,"/groups/irset/archives/SingleCell/projects/HumanRCC/ITH_ccRCC/201904/STEP4_DownstreamAnalysis_20190716/20190716_QCs_DFNorm/Aggr_ALL_ITH_ccRCC_sce_QCCellGenes_DF_20190718.rds")
##Load RDS:
#sce_QcCellsGenes_DF <- readRDS(file="/groups/irset/archives/SingleCell/projects/HumanRCC/ITH_ccRCC/201904/STEP4_DownstreamAnalysis_20190716/20190716_QCs_DFNorm/Aggr_ALL_ITH_ccRCC_sce_QCCellGenes_DF_20190718.rds")
####################################################################################

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])