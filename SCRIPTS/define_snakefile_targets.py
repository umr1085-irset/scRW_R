def dst(config, OUTDIR):
	'''
	Define Snakefile targets (completion files) from config file
	Parameters
	----------
	config : json
		Config parameters as a JSON object, generated from the config file
	OUTDIR : str
		Path of the outpout directory
	'''
	targets = [] # create empty list of target rules
	norms = []
	normalizations = config['NORMALIZATIONS'] # retrive all names of available normalizations
	if not any([config[x] for x in normalizations]): # retrieve all normalization methods True / False values
		if config['SEURAT_PIPE']:
			targets.append(config['STEP_SEURAT_PIPE_%s' % config['DEFAULT_NORMALIZATION']])
			norms.append(config['DEFAULT_NORMALIZATION'].lower())
		else:
			targets.append(config['STEP_%s' % config['DEFAULT_NORMALIZATION']]) # if all normalizations are False, use default normalization method
	else: # if one method is set to True
		for method in normalizations: # loop through methods
			if config[method]: # if method set to True
				norms.append(method.lower())
				if config['SEURAT_PIPE']:
					targets.append(config['STEP_SEURAT_PIPE_%s' % method])
				else:
					targets.append(config['STEP_%s' % method]) # append step
	targets.append('report')
	targets = [OUTDIR+".completion/%s" % x for x in targets] # update list of target paths with output directory and hidden completion folder
	
	return targets, norms

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