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
library(scran)

###############
# Load RDS file
###############
sce_QcCells = readRDS(file=snakemake@input[['rds_sce_cells']])

#######################
# Cell phase assignment
#######################
print('Retrieving cycle markers')
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran")) #list of 3 elements : "G1", "S" and "G2M"
rownames(sce_QcCells) <- rowData(sce_QcCells)$ID #Ensembl gene IDs as rownames

print('Computing cell cycle assignments')
assignments <- cyclone(sce_QcCells, pairs=hs.pairs) # assign a cell phase to each cell
colData(sce_QcCells)["CellCyclePhase"] <-  assignments$phases # store cell cycle phases in column data

#Save RDS:
print('Saving cell phase assignments')
saveRDS(assignments,snakemake@output[["rds_cell_phase"]]) # save assignments to file
saveRDS(sce_QcCells,snakemake@input[['rds_sce_cells']]) # update sce_QcCells object

###################
# Plot phase scores
###################
pdf(snakemake@output[['cell_cycle_plot']]) # create PDF plot
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=".")
dev.off()

#sce <- sce[,assignments$phases=="G1"] # do not grab G1 cells only for now (developmental data)
###########################
# Plot colored phase scores
###########################
ncells <- nrow(colData(sce_QcCells)) # total number of cells
is.G1 <- colData(sce_QcCells)$Barcode[colData(sce_QcCells)$CellCyclePhase=="G1"] # G1 cell phase
is.G1 <- na.omit(is.G1)
is.G2M <- colData(sce_QcCells)$Barcode[colData(sce_QcCells)$CellCyclePhase=="G2M"] # G2M cell phase
is.G2M <- na.omit(is.G2M)
is.S <- colData(sce_QcCells)$Barcode[colData(sce_QcCells)$CellCyclePhase=="S"] # S cell phase
is.S <- na.omit(is.S)

ALLcells <- sce_QcCells$Barcode # retrieve all cell barcodes
CellsAssigned <- c(is.G1,is.G2M,is.S) # retrieve all called cell barcodes
Cells_NoAssignment <- setdiff(ALLcells,CellsAssigned) # find cell barcodes with no assignment

cells.col <- rep("white",ncells) # default color
cells.col[ALLcells %in% as.vector(is.G1)] <- "red" # all G1 cells are assigned red
cells.col[ALLcells %in% as.vector(is.G2M)] <- "blue" # all G2M cells are assigned blue
cells.col[ALLcells %in% as.vector(is.S)] <- "green" # all S cells are assigned green
cells.col[ALLcells %in% Cells_NoAssignment] <- "black" # all unknown cycle cells are black

# define add_legend function for later use
add_legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

pdf(snakemake@output[['cell_cycle_plot_colored']]) # create PDF plot
plot(assignments$score$G1, assignments$score$G2M, col=cells.col, pch=".", xlab="G1 score", ylab="G2/M score")
add_legend("topright",legend=c("G1","G2M","S"), pch=20, col=c("red","blue","green"), horiz=TRUE, bty='n')
dev.off()

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])