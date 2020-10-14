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
library(DropletUtils)
library(scater)
library(SingleCellExperiment)

##################
# Define variables
################## 
OUTDIR = snakemake@params[["outdir"]] # grab output directory
PKGDIR = paste(normalizePath('.'),'/',sep='')
file.ribo.genes = paste(PKGDIR,'RSC/GENELISTS/',snakemake@params[["ribogenesfile"]],sep='')

######################
# Load aggregated data
######################
# Load H5 file
print('Loading aggregated data')
sce <- read10xCounts(snakemake@input[['aggrmatrix']], 
          sample.names = NULL,
          col.names = TRUE, type = "HDF5",
          version = "3", genome = NULL)

# remove dot and digit in gene IDs, rowData(sce)
g<-rowData(sce)$ID # grab genes names
g<-lapply(g, function(x) strsplit(x,'.',fixed = TRUE)[[1]][1]) # remove dot and digit based on split
rowData(sce)$ID <- g # replace ID column of rowData dataframe
row.names(sce) <- rowData(sce)$Symbol # replace row names with symbol

# correct Sample column in metadata
for(i in 1:length(snakemake@params[['samplelist']])){ # for index i in range of sample list
	sample = snakemake@params[['samplelist']][i] # get ith sample
	colData(sce)[grep(paste0('-',i), rownames(colData(sce))),]$Sample = sample # assign Sample column value to matching barcodes
}

# correct Sample column in metadata
#print(paste0('Running DoubletFinder on ', snakemake@params[['current_sample']])) # display current sample
#bcsgz=gzfile(snakemake@input[['bcsfile']],'rt') # barcodes file handle
#bcs=read.csv(bcsgz,header=F) # read barcodes
#bcs=bcs$V1 # extract barcodes to vector
#idx=match(snakemake@wildcards[["sample"]],snakemake@params[['samplelist']]) # get index of sample in aggregation sample list
#bcs=gsub("1", paste0(idx), bcs) # replace default "1" index in individual barcode list with sample index from aggregation list

########################
# Remove redundant genes
########################
print('Removing redundant genes')
sce <- sce[-which(duplicated(rowData(sce)$Symbol))] # filter out based on indices of rows with existing symbols

#########
# Cell QC
#########
if(snakemake@params[["species"]]=='HUMAN'){
	is.mito <- grep("^MT-",rownames(sce)) # find mitochondrial genes
} else {
	is.mito <- grep("^mt-",rownames(sce)) # find mitochondrial genes
}
ribo.data = read.table(file.ribo.genes, h=T) # read ribosomal gene list
is.ribo = which(rowData(sce)$ID %in% ribo.data$Gene_stable_ID) # find ribosomal genes

print('Computing cell QC metrics')
per.cell.qcmetrics = perCellQCMetrics(sce, subset=list(Mt=is.mito, Ri=is.ribo)) # compute per cell QC metrics

print('Computing feature QC metrics')
per.feature.qcmetrics = perFeatureQCMetrics(sce) # compute per feature QC metrics

colData(sce) <- cbind(colData(sce), per.cell.qcmetrics) # append per cell qc metrics to SCE object
rowData(sce) <- cbind(rowData(sce), per.feature.qcmetrics) # append per feature qc metrics to SCE object

##########
# Save RDS
##########
saveRDS(sce,snakemake@output[["rds_sce"]])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])