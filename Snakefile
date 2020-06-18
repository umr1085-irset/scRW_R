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

######################################
# CONFIG FILE
######################################
configfile: "config.yaml" # load config file

######################################
# PATHS
######################################
OUTDIR=os.path.join(config['outdir'],'') # grab folder information
INDIVDIR=os.path.join(config['indivdir'],'')

######################################
# TARGET
######################################
grp = [
	OUTDIR+".completion/step1", # create sce object
	OUTDIR+".completion/step2", # remove cell outliers
	OUTDIR+".completion/step3", # cell phase assignment
	OUTDIR+".completion/step4", # gene filtering
	OUTDIR+".completion/step5", # run DoubletFinder and filter out doublets
	#OUTDIR+".completion/step6", # deconvolution normalization
	OUTDIR+".completion/step6_merge" # individual deconvolution normalization + merge
]
rule all:
	input: grp

######################################
# SAMPLE LIST EXTRACTION
######################################
SAMPLES = sle(config['aggrfile'], config['sample_extraction_fromcol'])
check_samples(SAMPLES, config['indivdir'])

######################################
# PREPROCESSING
######################################
rule step1_create_sce_obj:
	input:
		aggrmatrix=config['aggrmatrix']
	output:
		rds_sce=OUTDIR+"objects/sce/sce.rds",
		step_complete=OUTDIR+".completion/step1"
	params:
		outdir=OUTDIR,
		ribogenesfile=config['ribogenesfile']
	script:
		"SCRIPTS/step1_create_sce_obj.R"

rule step2_cell_outliers:
	input:
		rds_sce=OUTDIR+"objects/sce/sce.rds"
	output:
		rds_sce_cells=OUTDIR+"objects/sce/sce_cells.rds",
		outliers_plot=OUTDIR+"QC/cell_outliers/QC_outlier_cells.pdf",
		step_complete=OUTDIR+".completion/step2"
	script:
		"SCRIPTS/step2_cell_outliers.R"

rule step3_cell_phase_assignment:
	input:
		rds_sce_cells=OUTDIR+"objects/sce/sce_cells.rds"
	output:
		rds_cell_phase=OUTDIR+"objects/cell_phase/cell_phase_assignments.rds",
		cell_cycle_plot=OUTDIR+"QC/cell_cyle/QC_CellCycleAssignment.pdf",
		cell_cycle_plot_colored=OUTDIR+"QC/cell_cyle/QC_CellCycleAssignment_colored.pdf",
		step_complete=OUTDIR+".completion/step3"
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
		QC_genes_AveGeneExpr_PCAOutliersRemoved_breaks100=OUTDIR+'QC/genes/QC_genes_AveGene_ExprPCAOutliersRemoved_breaks100.pdf.pdf',
		QC_genes_AveGeneExpr_PCAOutliersRemoved_LowAbundanceGenesRemoved_breaks100=OUTDIR+'QC/genes/QC_genes_AveGeneExpr_QcCellsGenes_breaks100.pdf',
		QC_genes_Top50Expr_GreyHist_QcCellsGenes=OUTDIR+'QC/genes/QC_genes_Top50Expr_GreyHist_QcCellsGenes.pdf',
		rds_sce_cells_genes=OUTDIR+"objects/sce/sce_cells_genes.rds",
		step_complete=OUTDIR+".completion/step4"
	script:
		"SCRIPTS/step4_gene_filtering.R"

rule step5_doublet_detection:
	input:
		rds_sce_cells_genes=OUTDIR+"objects/sce/sce_cells_genes.rds",
		bcsfile=INDIVDIR+"{sample}/filtered_feature_bc_matrix/barcodes.tsv.gz"
	params:
		samplelist=SAMPLES
	output:
		doublets_sample_file=OUTDIR+'DoubletFinder/{sample}_barcodes_doublets.txt',
		singlets_sample_file=OUTDIR+'DoubletFinder/{sample}_barcodes_singlets.txt',
		step_complete=OUTDIR+".completion/doubletfinder/{sample}"
	script:
		"SCRIPTS/step5_doublet_detection.R"

rule step5_doublet_filter:
	input:
		bcs_file_list=expand(OUTDIR+"DoubletFinder/{sample}_barcodes_singlets.txt",sample=SAMPLES),
		rds_sce_cells_genes=OUTDIR+"objects/sce/sce_cells_genes.rds",
	output:
		rds_sce_cells_genes_filtered=OUTDIR+"objects/sce/sce_cells_genes_filtered.rds",
		rds_sce_cells_genes_DF=OUTDIR+"objects/sce/sce_cells_genes_DF.rds",
		all_singlets=OUTDIR+'DoubletFinder/ALL_barcodes_singlets.txt',
		step_complete=OUTDIR+".completion/step5"
	script:
		"SCRIPTS/step5_doublet_filter.R"

rule step6_deconvolution_norm:
	input:
		rds_sce_cells_genes_filtered=OUTDIR+"objects/sce/sce_cells_genes_filtered.rds",
	output:
		rds_clusters=OUTDIR+"objects/clusters/ALL_clusters.rds",
		plot_Norm_HistSizeFactors=OUTDIR+"Normalization/Norm_HistSizeFactors.pdf",
		plot_Norm_SizeFactorsVsTotalCountsPerMillion=OUTDIR+"Normalization/Norm_SizeFactorsVsTotalCountsPerMillion.pdf",
		plot_Norm_SizeFactorsVsTotalCounts_smooth=OUTDIR+"Normalization/Norm_SizeFactorsVsTotalCounts_smooth.pdf",
		rds_sce_cells_genes_filtered_normed=OUTDIR+"objects/sce/sce_cells_genes_filtered_normed.rds",
		step_complete=OUTDIR+".completion/step6"
	params:
		individual_samples=0
	script:
		"SCRIPTS/step6_deconvolution_norm.R"

rule step6_individual_deconvolution_norm:
	input:
		sample_singlets=OUTDIR+"DoubletFinder/{sample}_barcodes_singlets.txt",
		rds_sce_cells_genes_filtered=OUTDIR+"objects/sce/sce_cells_genes_filtered.rds",
	output:
		rds_clusters=OUTDIR+"objects/clusters/individual_samples/{sample}_clusters.rds",
		plot_Norm_HistSizeFactors=OUTDIR+"Normalization/Indivual_samples/{sample}/Norm_HistSizeFactors.pdf",
		plot_Norm_SizeFactorsVsTotalCountsPerMillion=OUTDIR+"Normalization/Indivual_samples/{sample}/Norm_SizeFactorsVsTotalCountsPerMillion.pdf",
		plot_Norm_SizeFactorsVsTotalCounts_smooth=OUTDIR+"Normalization/Indivual_samples/{sample}/Norm_SizeFactorsVsTotalCounts_smooth.pdf",
		rds_sce_cells_genes_filtered_normed=OUTDIR+"objects/sce/individual_samples/{sample}/sce_cells_genes_filtered_normed.rds",
		step_complete=OUTDIR+".completion/individual_normalization/{sample}"
	params:
		individual_samples=1
	script:
		"SCRIPTS/step6_deconvolution_norm.R"

rule step6_merge_individual_normalizations:
	input:
		individual_files=expand(OUTDIR+"objects/sce/individual_samples/{sample}/sce_cells_genes_filtered_normed.rds",sample=SAMPLES)
	output:
		rds_sce_cells_genes_filtered_indivnormed_merged=OUTDIR+"objects/sce/sce_cells_genes_filtered_indivnormed_merged.rds",
		step_complete=OUTDIR+".completion/step6_merge"
	script:
		"SCRIPTS/step6_merge_individual_normalizations.R"

#rule step6_individual_deconvolution_norm:
#	input:
#		OUTDIR+"DoubletFinder/{sample}_barcodes_singlets.txt",
#		rds_sce_cells_genes_filtered=OUTDIR+"objects/sce/sce_cells_genes_filtered.rds",
#	output:
#		rds_clusters=OUTDIR+"objects/clusters/clusters.rds",
#		plot_Norm_HistSizeFactors=OUTDIR+"Normalization/Norm_HistSizeFactors.pdf",
#		plot_Norm_SizeFactorsVsTotalCountsPerMillion=OUTDIR+"Normalization/Norm_SizeFactorsVsTotalCountsPerMillion.pdf",
#		plot_Norm_SizeFactorsVsTotalCounts_smooth=OUTDIR+"Normalization/Norm_SizeFactorsVsTotalCounts_smooth.pdf",
#		rds_sce_cells_genes_filtered_normed=OUTDIR+"objects/sce/sce_cells_genes_filtered_normed.rds",
#		step_complete=OUTDIR+".completion/step6indidual"
#	script:
#		"SCRIPTS/step6_individual_deconvolution_norm.R"


