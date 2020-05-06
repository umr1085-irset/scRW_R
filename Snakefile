#############################################################
# Snakefile for single cell analysis workflow
# P. RIVAUD
# 2020/04
#############################################################

import os
import re
import sys
sys.path.append("SCRIPTS")
from sample_list_extraction import sle

######################################
# CONFIG FILE
######################################
configfile: "config.yaml" # load config file

######################################
# PATHS
######################################
OUTDIR=os.path.join(config['outdir'],'') # grab folder information
AGGRDIR=os.path.join(config['aggrdir'],'') # grab aggr folder info
AGGRMATRIX=config['aggrmatrix'] # grab aggr folder info
AGGRFILE=config['aggrfile'] # grab aggr folder info
RIBOGENESFILE=config['ribogenesfile'] # grab aggr folder info

######################################
# TARGET
######################################
grp = [
		OUTDIR+".completion/step1",
		OUTDIR+".completion/step2",
		OUTDIR+".completion/step3",
		OUTDIR+".completion/step4"
	]
rule all:
	input: grp

######################################
# SAMPLE LIST EXTRACTION
######################################
SAMPLES = sle(AGGRFILE)
print(SAMPLES)

######################################
# PREPROCESSING
######################################
rule step1_create_sce_obj:
	input:
		aggrfile=AGGRMATRIX
	output:
		rds_sce=OUTDIR+"objects/sce/sce.rds",
		step_complete=OUTDIR+".completion/step1"
	params:
		outdir=OUTDIR,
		ribogenesfile=RIBOGENESFILE
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
		#rds_clusters=OUTDIR+"objects/clusters/clusters.rds",
		step_complete=OUTDIR+".completion/step4"
	script:
		"SCRIPTS/step4_gene_filtering.R"

'''
rule step5_doublet_detection:
	input:
'''