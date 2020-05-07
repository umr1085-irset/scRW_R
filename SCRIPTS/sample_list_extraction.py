#############################################################
# Script for single cell analysis workflow
# P. RIVAUD
# 2020/04
#############################################################

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
		return df.library_id.values.tolist()
	elif fromcol=='molecule_h5':
		return [x.split('/')[-3] for x in df.molecule_h5]
	else:
		raise ValueError('Invalid value for fromcol parameter. Must be one of library_id, molecule_h5')

