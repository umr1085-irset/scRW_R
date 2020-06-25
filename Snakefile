#############################################################
# Snakefile for single cell analysis workflow
# P. RIVAUD
# 2020/04
#############################################################

import os
import re
import sys
sys.path.append("SCRIPTS")
from sample_list_extraction import sle, check_samples

from snakemake.logging import logger

######################################
# CONFIG FILE
######################################
configfile: "config.yaml" # load config file

######################################
# PATHS
######################################
OUTDIR=os.path.join(config['OUTDIR'],'') # grab folder information
INDIVDIR=os.path.join(config['INDIVDIR'],'')

######################################
# TARGETS
######################################
# steps that should be computed automatically:
# - step1_create_sce_obj
# - step2_cell_outliers
# - step4_gene_filtering
# - step5_doublet_filter

targets = [] # create empty list of target rules
for step in config['OPTIONAL_STEPS']: # loop through optional steps
	if config[step]: targets.append(config['STEP_%s' % step]) # if step is True, add step completion file to target list

normalizations = config['NORMALIZATIONS'] # retrive all names of available normalizations
if not any([config[x] for x in normalizations]): # retrieve all normalization methods True / False values
	targets.append(config['STEP_%s' % config['DEFAULT_NORMALIZATION']]) # if all normalizations are False, use default normalization method
else: # if one method is set to True
	for method in normalizations: # loop through methods
		if config[method]: # if method set to True
			targets.append(config['STEP_%s' % method]) # append step

targets = [OUTDIR+".completion/%s" % x for x in targets] # update list of target paths with output directory and hidden completion folder

#print(targets)
rule all:
	input: targets

######################################
# SAMPLE LIST EXTRACTION
######################################
SAMPLES = sle(config['AGGRFILE'], config['SAMPLE_EXTRACTION_FROMCOL'])
check_samples(SAMPLES, config['INDIVDIR'])

######################################
# PREPROCESSING
######################################
rule step1_create_sce_obj:
	input:
		aggrmatrix=config['AGGRMATRIX']
	output:
		rds_sce=OUTDIR+"objects/sce/sce.rds",
		step_complete=OUTDIR+".completion/step1_create_sce_obj"
	params:
		outdir=OUTDIR,
		ribogenesfile=config['RIBOGENESFILE']
	script:
		"SCRIPTS/step1_create_sce_obj.R"

rule step2_cell_outliers:
	input:
		rds_sce=OUTDIR+"objects/sce/sce.rds"
	output:
		rds_sce_cells=OUTDIR+"objects/sce/sce_cells.rds",
		outliers_plot=OUTDIR+"QC/cell_outliers/QC_outlier_cells.pdf",
		step_complete=OUTDIR+".completion/step2_cell_outliers"
	script:
		"SCRIPTS/step2_cell_outliers.R"

rule step3_cell_phase_assignment:
	input:
		rds_sce_cells=OUTDIR+"objects/sce/sce_cells.rds"
	output:
		rds_cell_phase=OUTDIR+"objects/cell_phase/cell_phase_assignments.rds",
		cell_cycle_plot=OUTDIR+"QC/cell_cyle/QC_CellCycleAssignment.pdf",
		cell_cycle_plot_colored=OUTDIR+"QC/cell_cyle/QC_CellCycleAssignment_colored.pdf",
		step_complete=OUTDIR+".completion/step3_cell_phase_assignment"
	script:
		"SCRIPTS/step3_cell_phase_assignment.R"

rule step4_gene_filtering:
	input:
		rds_sce=OUTDIR+"objects/sce/sce.rds",
		rds_sce_cells=OUTDIR+"objects/sce/sce_cells.rds"
	output:
		QC_genes_AveCounts_PCAOutliersRemoved=OUTDIR+"QC/genes/QC_genes_AveCounts_PCAOutliersRemoved.pdf",
		QC_genes_NumCells_Counts_PCAOutliersRemoved=OUTDIR+"QC/genes/QC_genes_NumCells_Counts_PCAOutliersRemoved.pdf",
		QC_genes_AveGeneExpr_breaks100=OUTDIR+'QC/genes/QC_genes_AveGeneExpr_breaks100.pdf',
		QC_genes_AveGeneExpr_PCAOutliersRemoved_breaks100=OUTDIR+'QC/genes/QC_genes_AveGeneExpr_PCAOutliersRemoved_breaks100.pdf',
		QC_genes_AveGeneExpr_PCAOutliersRemoved_LowAbundanceGenesRemoved_breaks100=OUTDIR+'QC/genes/QC_genes_AveGeneExpr_QcCellsGenes_breaks100.pdf',
		QC_genes_Top50Expr_GreyHist_QcCellsGenes=OUTDIR+'QC/genes/QC_genes_Top50Expr_GreyHist_QcCellsGenes.pdf',
		rds_sce_cells_genes=OUTDIR+"objects/sce/sce_cells_genes.rds",
		step_complete=OUTDIR+".completion/step4_gene_filtering"
	script:
		"SCRIPTS/step4_gene_filtering.R"

rule step5_doublet_detection:
	input:
		rds_sce_cells_genes=OUTDIR+"objects/sce/sce_cells_genes.rds",
		bcsfile=INDIVDIR+"{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
	params:
		samplelist=SAMPLES,
		current_sample="{sample}"
	output:
		doublets_sample_file=OUTDIR+'DoubletFinder/barcode_lists/{sample}_barcodes_doublets.txt',
		singlets_sample_file=OUTDIR+'DoubletFinder/barcode_lists/{sample}_barcodes_singlets.txt',
		step_complete=OUTDIR+".completion/doubletfinder/{sample}"
	script:
		"SCRIPTS/step5_doublet_detection.R"

rule step5_doublet_filter:
	input:
		bcs_file_list=expand(OUTDIR+"DoubletFinder/barcode_lists/{sample}_barcodes_singlets.txt",sample=SAMPLES),
		rds_sce_cells_genes=OUTDIR+"objects/sce/sce_cells_genes.rds",
	output:
		rds_sce_cells_genes_singlets=OUTDIR+"objects/sce/sce_cells_genes_singlets.rds",
		rds_sce_cells_genes_DF=OUTDIR+"objects/sce/sce_cells_genes_DF.rds",
		all_singlets=OUTDIR+'DoubletFinder/barcode_lists/ALL_barcodes_singlets.txt',
		step_complete=OUTDIR+".completion/step5_doublet_filter"
	script:
		"SCRIPTS/step5_doublet_filter.R"

rule step6_norm_scran_deconvolution:
	input:
		rds_sce_cells_genes_singlets=OUTDIR+"objects/sce/sce_cells_genes_singlets.rds",
	output:
		rds_clusters=OUTDIR+"Normalization/scran_deconvolution/ALL_normalization_clusters.rds",
		plot_Norm_HistSizeFactors=OUTDIR+"Normalization/scran_deconvolution/Norm_HistSizeFactors.pdf",
		plot_Norm_SizeFactorsVsTotalCountsPerMillion=OUTDIR+"Normalization/scran_deconvolution/Norm_SizeFactorsVsTotalCountsPerMillion.pdf",
		plot_Norm_SizeFactorsVsTotalCounts_smooth=OUTDIR+"Normalization/scran_deconvolution/Norm_SizeFactorsVsTotalCounts_smooth.pdf",
		rds_sce_cells_genes_singlets_normed=OUTDIR+"objects/sce/sce_cells_genes_singlets_scran_deconvolution.rds",
		step_complete=OUTDIR+".completion/step6_scran_deconvolution"
	params:
		individual_samples=0
	script:
		"SCRIPTS/step6_norm_scran_deconvolution.R"

rule step6_norm_scran_deconvolution_individual:
	input:
		sample_singlets=OUTDIR+"DoubletFinder/barcode_lists/{sample}_barcodes_singlets.txt",
		rds_sce_cells_genes_singlets=OUTDIR+"objects/sce/sce_cells_genes_singlets.rds",
	output:
		rds_clusters=OUTDIR+"Normalization/scran_deconvolution/individual_samples/{sample}/normalization_clusters.rds",
		plot_Norm_HistSizeFactors=OUTDIR+"Normalization/scran_deconvolution/individual_samples/{sample}/Norm_HistSizeFactors.pdf",
		plot_Norm_SizeFactorsVsTotalCountsPerMillion=OUTDIR+"Normalization/scran_deconvolution/individual_samples/{sample}/Norm_SizeFactorsVsTotalCountsPerMillion.pdf",
		plot_Norm_SizeFactorsVsTotalCounts_smooth=OUTDIR+"Normalization/scran_deconvolution/individual_samples/{sample}/Norm_SizeFactorsVsTotalCounts_smooth.pdf",
		rds_sce_cells_genes_singlets_normed=OUTDIR+"objects/sce/individual_samples/{sample}/sce_cells_genes_singlets_scran_deconvolution.rds",
		step_complete=OUTDIR+".completion/individual_normalization/step6_scran_deconvolution_{sample}"
	params:
		individual_samples=1
	script:
		"SCRIPTS/step6_norm_scran_deconvolution.R"

rule step6_merge_individual_scran_deconvolutions:
	input:
		individual_files=expand(OUTDIR+"objects/sce/individual_samples/{sample}/sce_cells_genes_singlets_scran_deconvolution.rds", sample=SAMPLES)
	output:
		rds_sce_cells_genes_singlets_scran_deconvolution_merged=OUTDIR+"objects/sce/sce_cells_genes_singlets_scran_deconvolution_merged.rds",
		step_complete=OUTDIR+".completion/step6_scran_deconvolution_merged"
	script:
		"SCRIPTS/step6_merge_individual_scran_deconvolutions.R"

rule step6_seurat_sctransform:
	input:
		rds_sce_cells_genes_singlets=OUTDIR+"objects/sce/sce_cells_genes_singlets.rds"
	output:
		rds_sce_cells_genes_singlets_normed=OUTDIR+"objects/sce/sce_cells_genes_singlets_seurat_sctransform.rds",
		step_complete=OUTDIR+".completion/step6_seurat_sctransform"
	script:
		"SCRIPTS/step6_seurat_sctransform.R"



###################################################################################################################################################
###################################################################################################################################################
###################################################################################################################################################
'''
targets = ['step1_create_sce_obj', 'step2_cell_outliers']
if config['CELL_PHASE_ASSIGNMENT']: targets.append('step3_cell_phase_assignment')
targets += ['step4_gene_filtering', 'step5_doublet_filter']
normalizations = config['NORMALIZATIONS']
normalizations = [config[x] for x in normalizations]
if not any(normalizations): default = True
if config['SCRAN_DECONVOLUTION']
if config['SCRAN_DECONVOLUTION_INDIVIDUAL']
targets = [OUTDIR+".completion/%s" % x for x in targets]
'''
'''
targets = [
	#OUTDIR+".completion/step1_create_sce_obj", # create sce object
	#OUTDIR+".completion/step2_cell_outliers", # remove cell outliers
	#OUTDIR+".completion/step3_cell_phase_assignment", # cell phase assignment - optional
	#OUTDIR+".completion/step4_gene_filtering", # gene filtering
	OUTDIR+".completion/step5_doublet_filter", # run DoubletFinder and filter out doublets
	#OUTDIR+".completion/step6_scran_deconvolution", # deconvolution normalization
	#OUTDIR+".completion/step6_scran_deconvolution_merged", # individual deconvolution normalization + merge
]
'''