#############################################################
# Script for single cell analysis workflow
# P. RIVAUD
# 2020/04
#############################################################

import os
import pandas as pd

def sle(filename,fromcol='library_id'):
	'''
	Extract sample list from aggregation file used in cellranger aggr pipeline
	Parameters
	----------
	fromcol : str
		Column from which to extract the sample names. Can be `library_id` or `molecule_h5`
	'''
	df = pd.read_csv(filename, header=0)
	if fromcol=='library_id':
		return df.library_id.values.tolist() # sample names are the library ids
	elif fromcol=='molecule_h5':
		return [x.split('/')[-3] for x in df.molecule_h5] # sample names are the folder names in the aggregation step molecule.h5 paths
	else:
		raise ValueError('Invalid value for fromcol parameter. Must be one of library_id, molecule_h5')

def check_samples(L,D):
	'''
	Check if samples in `L`exist in `D`
	Parameters
	----------
	L : list
		List of sample names
	D : str
		Directory of individual (unaggregated) samples
	'''
	dirlist = next(os.walk(D))[1] # get list of samples (folders) in unaggregated data folder
	if set(L).issubset(dirlist) == False:
		raise Exception('Missing some individual samples in %s folder' % D)