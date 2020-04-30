# Load matrix from folder
#library(Seurat)
print('Loading matrix')
#cellbarcodes <- read.table("/Users/paulrivaud/Documents/data/20200306/SINS/SIN1/filtered_feature_bc_matrix/barcodes.tsv.gz")
#genenames    <- read.table("/Users/paulrivaud/Documents/data/20200306/SINS/SIN1/filtered_feature_bc_matrix/features_modified.tsv.gz")
#molecules    <- Matrix::readMM("/Users/paulrivaud/Documents/data/20200306/SINS/SIN1/filtered_feature_bc_matrix/matrix.mtx.gz")

# work with Ensembl gene IDs genenames$V1, OR genenames$V2 if want to work with geneNames
rownames(molecules) <- genenames$V2
colnames(molecules) <- cellbarcodes$V1
print('matrix dim:')
print(dim(molecules))


# Build a scater SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData=colnames(molecules),rowData=rownames(molecules))
sce2 <- read10xCounts(samples, sample.names = names(samples),
  col.names = TRUE, type = "HDF5",
  version = "3", genome = NULL)

# Expression data in sce object:
print('sce obj dim:')
print(dim(counts(sce)))
colData(sce)
colData(sce,internal=TRUE)
#DataFrame with 147060 rows and 1 column --> cells
### Add GeneIDs and GeneNames informations in the sce object ###
rowData(sce)["GeneID"]   <- genenames$V1
rowData(sce)["GeneName"] <- genenames$V2

assays(sce)
#List of length 1
#names(1): counts
assayNames(sce) #"counts"

colData(sce)
sce$X
# -> Liste des n cellules
# sce$name avec name = name of a column of colData