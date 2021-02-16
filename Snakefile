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
from define_snakefile_targets import dst

######################################
# CONFIG FILE
######################################
configfile: "config.yaml" # load config file

######################################
# METADATA
######################################
if config['SPECIES'].upper()=='HUMAN':
	ribogenesfile=config['HUMAN_RIBOGENESFILE']
elif config['SPECIES'].upper()=='MOUSE':
	ribogenesfile=config['MOUSE_RIBOGENESFILE']
else:
	raise ValueError('SPECIES in config file must be HUMAN or MOUSE. Exiting')
regressoncellcyles=1 if config['REGRESSCELLCYCLES'] else 0

######################################
# PATHS
######################################
OUTDIR=os.path.join(os.path.abspath(config['OUTDIR']),'') # grab folder information
INDIVDIR=os.path.join(config['INDIVDIR'],'')
RSCDIR=os.path.join(os.getcwd(),'RSC','')

######################################
# TARGETS
######################################
target_all = OUTDIR+".completion/report" if config['REPORT'] else OUTDIR+".completion/gather"
rule all:
	input: target_all

targets, norms = dst(config, OUTDIR) # Define Snakefile Targets
rule gather:
	input: targets
	output: OUTDIR+".completion/gather"

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
		datalog=OUTDIR+"report/datalog/datalogstep1_num_genes_removed.txt",
		datalog_df=OUTDIR+"report/datalog/datalogstep1_df.csv",
		step_complete=OUTDIR+".completion/step1_create_sce_obj"
	params:
		ribogenesfile=RSCDIR+'GENELISTS/%s' % ribogenesfile,
		species=config['SPECIES'].upper(),
		samplelist=SAMPLES
	script:
		"SCRIPTS/step1_create_sce_obj.R"

rule step2_cell_outliers:
	input:
		rds_sce=OUTDIR+"objects/sce/sce.rds"
	output:
		rds_sce_cells=OUTDIR+"objects/sce/sce_cells.rds",
		outliers_plot=OUTDIR+"QC/cell_outliers/QC_outlier_cells.pdf",
		datalog_df=OUTDIR+"report/datalog/datalogstep2_df.csv",
		step_complete=OUTDIR+".completion/step2_cell_outliers"
	script:
		"SCRIPTS/step2_cell_outliers.R"

rule step3_cell_phase_assignment:
	input:
		rds_sce_cells=OUTDIR+"objects/sce/sce_cells.rds"
	output:
		rds_cell_phase=OUTDIR+"objects/sce/sce_cells_cellphase.rds",
		cell_cycle_plot=OUTDIR+"QC/cell_cycle/QC_CellCycleAssignment.pdf",
		cell_cycle_plot_colored=OUTDIR+"QC/cell_cycle/QC_CellCycleAssignment_colored.pdf",
		step_complete=OUTDIR+".completion/step3_cell_phase_assignment"
	params:
		species=config['SPECIES'].upper()
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
		datalog_df=OUTDIR+"report/datalog/datalogstep4_df.csv",
		step_complete=OUTDIR+".completion/step4_gene_filtering"
	script:
		"SCRIPTS/step4_gene_filtering.R"

rule step5_doublet_detection:
	input:
		rds_sce_cells_genes=OUTDIR+"objects/sce/sce_cells_genes.rds",
		bcsfile=os.path.join(INDIVDIR+"{sample}",config['SUBPATH_IN_INDIVDIRS'],"barcodes.tsv.gz")
		#bcsfile=INDIVDIR+"{sample}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
	params:
		samplelist=SAMPLES,
		current_sample="{sample}"
	output:
		doublets_sample_file=OUTDIR+'DoubletFinder/barcode_lists/{sample}_barcodes_doublets.txt',
		singlets_sample_file=OUTDIR+'DoubletFinder/barcode_lists/{sample}_barcodes_singlets.txt',
		step_complete=OUTDIR+".completion/doubletfinder/{sample}"
	resources:
		load=1
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
		rds_clusters=OUTDIR+"normalization/scran_deconvolution/ALL_normalization_clusters.rds",
		plot_Norm_HistSizeFactors=OUTDIR+"normalization/scran_deconvolution/Norm_HistSizeFactors.pdf",
		plot_Norm_SizeFactorsVsTotalCountsPerMillion=OUTDIR+"normalization/scran_deconvolution/Norm_SizeFactorsVsTotalCountsPerMillion.pdf",
		plot_Norm_SizeFactorsVsTotalCounts_smooth=OUTDIR+"normalization/scran_deconvolution/Norm_SizeFactorsVsTotalCounts_smooth.pdf",
		rds_sce_cells_genes_singlets_normed=OUTDIR+"objects/sce/sce_cells_genes_singlets_scran_deconvolution.rds",
		step_complete=OUTDIR+".completion/step6_norm_scran_deconvolution"
	params:
		individual_samples=0
	script:
		"SCRIPTS/step6_norm_scran_deconvolution.R"

rule step6_norm_seurat_sctransform:
	input:
		rds_sce_cells_genes_singlets=OUTDIR+"objects/sce/sce_cells_genes_singlets.rds"
	output:
		seurat_cells_genes_singlets_normed=OUTDIR+"objects/seurat/seurat_cells_genes_singlets_seurat_sctransform.rds",
		step_complete=OUTDIR+".completion/step6_norm_seurat_sctransform"
	params:
		regressoncellcyles=regressoncellcyles,
		future_globals_maxsize=config['FUTURE_GLOBALS_MAXSIZE']
	script:
		"SCRIPTS/step6_norm_seurat_sctransform.R"

rule step6_norm_seurat_lognorm:
	input:
		rds_sce_cells_genes_singlets=OUTDIR+"objects/sce/sce_cells_genes_singlets.rds"
	output:
		seurat_cells_genes_singlets_normed=OUTDIR+"objects/seurat/seurat_cells_genes_singlets_seurat_lognorm.rds",
		step_complete=OUTDIR+".completion/step6_norm_seurat_lognorm"
	script:
		"SCRIPTS/step6_norm_seurat_lognorm.R"

rule step7_seurat_pipe_scran_deconvolution:
	input:
		rds_sce_cells_genes_singlets_normed=OUTDIR+"objects/sce/sce_cells_genes_singlets_scran_deconvolution.rds"
	output:
		rds_seurat=OUTDIR+'objects/seurat/seurat_cells_genes_singlets_scran_deconvolution_dimred.rds',
		plot_pca_clusters=OUTDIR+'normalization/scran_deconvolution/pca_clusters_normed_scran_deconvolution.pdf',
		plot_pca_cellphase=OUTDIR+'normalization/scran_deconvolution/pca_cellphase_normed_scran_deconvolution.pdf',
		plot_umap_clusters=OUTDIR+'normalization/scran_deconvolution/umap_clusters_normed_scran_deconvolution.pdf',
		plot_umap_cellphase=OUTDIR+'normalization/scran_deconvolution/umap_cellphase_normed_scran_deconvolution.pdf',
		datalog_df=OUTDIR+'report/datalog/datalogstep7_scran_deconvolution_df.csv',
		step_complete=OUTDIR+".completion/step7_seurat_pipe_scran_deconvolution"
	params:
		seuratinput=0,
		scaled=0,
		use_jack_straw=1,
		regressoncellcyles=regressoncellcyles
	script:
		"SCRIPTS/step7_seurat3_pipe.R"

rule step7_seurat_pipe_seurat_sctransform:
	input:
		rds_sce_cells_genes_singlets_normed=OUTDIR+"objects/seurat/seurat_cells_genes_singlets_seurat_sctransform.rds"
	output:
		rds_seurat=OUTDIR+'objects/seurat/seurat_cells_genes_singlets_seurat_sctransform_dimred.rds',
		plot_pca_clusters=OUTDIR+'normalization/seurat_sctransform/pca_clusters_seurat_sctransform.pdf',
		plot_pca_cellphase=OUTDIR+'normalization/seurat_sctransform/pca_cellphase_normed_seurat_sctransform.pdf',
		plot_umap_clusters=OUTDIR+'normalization/seurat_sctransform/umap_clusters_normed_seurat_sctransform.pdf',
		plot_umap_cellphase=OUTDIR+'normalization/seurat_sctransform/umap_cellphase_normed_seurat_sctransform.pdf',
		datalog_df=OUTDIR+'report/datalog/datalogstep7_seurat_sctransform_df.csv',
		step_complete=OUTDIR+".completion/step7_seurat_pipe_seurat_sctransform"
	params:
		seuratinput=1,
		scaled=1,
		use_jack_straw=0,
		regressoncellcyles=regressoncellcyles
	script:
		"SCRIPTS/step7_seurat3_pipe.R"

rule step7_seurat_pipe_seurat_lognorm:
	input:
		rds_sce_cells_genes_singlets_normed=OUTDIR+"objects/seurat/seurat_cells_genes_singlets_seurat_lognorm.rds"
	output:
		rds_seurat=OUTDIR+'objects/seurat/seurat_cells_genes_singlets_seurat_lognorm_dimred.rds',
		plot_pca_clusters=OUTDIR+'normalization/seurat_lognorm/pca_clusters_normed_seurat_lognorm.pdf',
		plot_pca_cellphase=OUTDIR+'normalization/seurat_lognorm/pca_cellphase_normed_seurat_lognorm.pdf',
		plot_umap_clusters=OUTDIR+'normalization/seurat_lognorm/umap_clusters_normed_seurat_lognorm.pdf',
		plot_umap_cellphase=OUTDIR+'normalization/seurat_lognorm/umap_cellphase_normed_seurat_lognorm.pdf',
		datalog_df=OUTDIR+'report/datalog/datalogstep7_seurat_lognorm_df.csv',
		step_complete=OUTDIR+".completion/step7_seurat_pipe_seurat_lognorm",
	params:
		seuratinput=1,
		scaled=0,
		use_jack_straw=1,
		regressoncellcyles=regressoncellcyles
	script:
		"SCRIPTS/step7_seurat3_pipe.R"

rule create_report:
	input:
		plotlyjs_file=RSCDIR+"WEB/plotly-latest.min.js",
		datalog_step1_genes_rm = OUTDIR+"report/datalog/datalogstep1_num_genes_removed.txt",
		datalog_df_step1 = OUTDIR+"report/datalog/datalogstep1_df.csv",
		datalog_df_step2 = OUTDIR+"report/datalog/datalogstep2_df.csv",
		datalog_df_step4 = OUTDIR+"report/datalog/datalogstep4_df.csv",
		steps_wait=OUTDIR+".completion/gather"
	params:
		css_file='https://bootswatch.com/4/spacelab/bootstrap.min.css',
		samplelist=SAMPLES,
		html_template_dir=RSCDIR+"WEB/jinja_templates",
		species=config['SPECIES'].upper(),
		basedir_step7 = OUTDIR+'report/datalog/',
		nsamples = len(SAMPLES),
		samples=SAMPLES,
		project_name=config['PROJECT_NAME'],
		doubletfinderdir=OUTDIR+'DoubletFinder/barcode_lists'
	output:
		html_report = OUTDIR+'report/report.html',
		step_complete=OUTDIR+".completion/report"
	script:
		"SCRIPTS/create_report.py"

#rule step6_norm_scran_deconvolution_individual:
#	input:
#		sample_singlets=OUTDIR+"DoubletFinder/barcode_lists/{sample}_barcodes_singlets.txt",
#		rds_sce_cells_genes_singlets=OUTDIR+"objects/sce/sce_cells_genes_singlets.rds",
#	output:
#		rds_clusters=OUTDIR+"normalization/scran_deconvolution/individual_samples/{sample}/normalization_clusters.rds",
#		plot_Norm_HistSizeFactors=OUTDIR+"normalization/scran_deconvolution/individual_samples/{sample}/Norm_HistSizeFactors.pdf",
#		plot_Norm_SizeFactorsVsTotalCountsPerMillion=OUTDIR+"normalization/scran_deconvolution/individual_samples/{sample}/Norm_SizeFactorsVsTotalCountsPerMillion.pdf",
#		plot_Norm_SizeFactorsVsTotalCounts_smooth=OUTDIR+"normalization/scran_deconvolution/individual_samples/{sample}/Norm_SizeFactorsVsTotalCounts_smooth.pdf",
#		rds_sce_cells_genes_singlets_normed=OUTDIR+"objects/sce/individual_samples/{sample}/sce_cells_genes_singlets_scran_deconvolution.rds",
#		step_complete=OUTDIR+".completion/individual_normalization/step6_scran_deconvolution_{sample}"
#	params:
#		individual_samples=1
#	script:
#		"SCRIPTS/step6_norm_scran_deconvolution.R"
#
#rule step6_merge_individual_scran_deconvolutions:
#	input:
#		individual_files=expand(OUTDIR+"objects/sce/individual_samples/{sample}/sce_cells_genes_singlets_scran_deconvolution.rds", sample=SAMPLES)
#	output:
#		rds_sce_cells_genes_singlets_scran_deconvolution_merged=OUTDIR+"objects/sce/sce_cells_genes_singlets_scran_deconvolution_merged.rds",
#		step_complete=OUTDIR+".completion/step6_scran_deconvolution_merged"
#	script:
#		"SCRIPTS/step6_merge_individual_scran_deconvolutions.R"

#rule step7_seurat_pipe_scran_deconvolution_merged:
#	input:
#		rds_sce_cells_genes_singlets_normed=OUTDIR+"objects/sce/sce_cells_genes_singlets_scran_deconvolution_merged.rds"
#	output:
#		rds_seurat=OUTDIR+'objects/seurat/sce_cells_genes_singlets_scran_deconvolution_merged_dimred.rds',
#		step_complete=OUTDIR+".completion/step7_seurat_pipe_scran_deconvolution_merged"
#	params:
#		seuratinput=0
#	script:
#		"SCRIPTS/step7_seurat3_pipe.R"