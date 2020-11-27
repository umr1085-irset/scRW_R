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

########################
# Remove redundant genes
########################
print('Removing redundant genes')
number_genes_removed <- length(which(duplicated(rowData(sce)$Symbol)))
sce <- sce[-which(duplicated(rowData(sce)$Symbol))] # filter out based on indices of rows with existing symbols

#########
# Cell QC
#########
if(snakemake@params[["species"]]=='HUMAN'){
	is.mito <- grep("^MT-",rownames(sce)) # find mitochondrial genes
} else {
	is.mito <- grep("^mt-",rownames(sce)) # find mitochondrial genes
}
ribo.data = read.table(snakemake@params[["ribogenesfile"]], h=T) # read ribosomal gene list
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

##########
# Log data
##########
write(number_genes_removed, snakemake@output[["datalog"]])
write.csv(colData(sce)[c('sum','detected','Sample')], snakemake@output[["datalog_df"]])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])

#########################
#########################
## R tests
#########################
#########################
#library(DropletUtils)
#library(scater)
#library(SingleCellExperiment)
#
#sce <- read10xCounts('TESTDATA/small_aggr/filtered_feature_bc_matrix.h5', 
#          sample.names = NULL,
#          col.names = TRUE, type = "HDF5",
#          version = "3", genome = NULL)
#
#
##dataloglines <- c(paste0('total_num_cells,',dim(sce)[2]))
#json = list()
#json['total_num_cells'] = dim(sce)[2]
#
#g<-rowData(sce)$ID # grab genes names
#g<-lapply(g, function(x) strsplit(x,'.',fixed = TRUE)[[1]][1]) # remove dot and digit based on split
#rowData(sce)$ID <- g # replace ID column of rowData dataframe
#row.names(sce) <- rowData(sce)$Symbol # replace row names with symbol
#
##['m_10+2_SIN1', 'f_9+4_SIN2']
## correct Sample column in metadata
#for(i in 1:length(c('m_10+2_SIN1', 'f_9+4_SIN2'))){ # for index i in range of sample list
#	sample = c('m_10+2_SIN1', 'f_9+4_SIN2')[i] # get ith sample
#	colData(sce)[grep(paste0('-',i), rownames(colData(sce))),]$Sample = sample # assign Sample column value to matching barcodes
#}
#
#dataloglines <-c(dataloglines, paste0('num_genes_removed,',length(which(duplicated(rowData(sce)$Symbol)))))
#sce <- sce[-which(duplicated(rowData(sce)$Symbol))] # filter out based on indices of rows with existing symbols
#
#is.mito <- grep("^MT-",rownames(sce))
#ribo.data = read.table('RSC/GENELISTS/mart_export_RiboGenes_GO0005840_human.txt', h=T) # read ribosomal gene list
#is.ribo = which(rowData(sce)$ID %in% ribo.data$Gene_stable_ID) # find ribosomal genes
#
#print('Computing cell QC metrics')
#per.cell.qcmetrics = perCellQCMetrics(sce, subset=list(Mt=is.mito, Ri=is.ribo)) # compute per cell QC metrics
#
#print('Computing feature QC metrics')
#per.feature.qcmetrics = perFeatureQCMetrics(sce) # compute per feature QC metrics
#
#colData(sce) <- cbind(colData(sce), per.cell.qcmetrics) # append per cell qc metrics to SCE object
#rowData(sce) <- cbind(rowData(sce), per.feature.qcmetrics) # append per feature qc metrics to SCE object
#
#write.csv(colData(sce)[c('sum','detected','Sample')], "OUTPUT/report/datalog/datalogstep1_df.csv")
#
#json = list(samples=samples,meta=json)
#datalog <- file("OUTPUT/report/datalog/datalogstep1")
#writeLines(dataloglines, datalog)
#close(datalog)#