#############################################################
# PREPROCESSING SCRIPT
# P. RIVAUD
# 2020/04
#############################################################

# In R, Snakefile parameters can be accessed with the following:
# print(snakemake@input[[1]])
# print(snakemake@output[["rds_sce"]])
# file.create(snakemake@output[["rds_sce"]])

#print(snakemake@input[['individual_files']])
sce_list = lapply(snakemake@input[['individual_files']], function(x) readRDS(file=x))
sce_merged = do.call(cbind, sce_list)

#############
# Save to RDS
#############
saveRDS(sce_merged, snakemake@output[['rds_sce_cells_genes_singlets_scran_deconvolution_merged']])

###############
# Complete step
###############
file.create(snakemake@output[["step_complete"]])