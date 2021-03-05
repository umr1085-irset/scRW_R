import scipy.io as sio
import scipy.sparse as ss
import gzip
import os
import numpy as np
import shutil

def replace1(s,n):
    return s.replace('1',n)

def read_barcodes(path,i):
    f=gzip.open(path,'rt')
    bcs=f.read().strip()
    bcs = np.array(bcs.split('\n'))
    f.close()
    if i!=1:
        bcs = np.vectorize(replace1)(bcs,str(i))
    return bcs

def compute_aggregation(config,samples):
	OUTDIR=os.path.join(os.path.abspath(config['OUTDIR']),'')
	mtx_path = os.path.join(config['INDIVDIR'], samples[0], config['SUBPATH_IN_INDIVDIRS'], 'matrix.mtx.gz')
	M = sio.mmread(mtx_path).tocsc()
	BCS = read_barcodes(os.path.join(config['INDIVDIR'], samples[0], config['SUBPATH_IN_INDIVDIRS'], 'barcodes.tsv.gz'),1)

	for i,sample in enumerate(samples[1:]):
		idx = i+2
		mtx_path = os.path.join(config['INDIVDIR'], sample, config['SUBPATH_IN_INDIVDIRS'], 'matrix.mtx.gz')
		M = ss.hstack([M,sio.mmread(mtx_path).tocsc()])
		BCS = np.concatenate((BCS,read_barcodes(os.path.join(config['INDIVDIR'], sample, config['SUBPATH_IN_INDIVDIRS'], 'barcodes.tsv.gz'),idx)))

	os.makedirs(OUTDIR, exist_ok=True)
	os.makedirs(os.path.join(OUTDIR,'aggrdata'), exist_ok=True)

	# write matrix
	sio.mmwrite(os.path.join(OUTDIR,'aggrdata','matrix.mtx'),M)
	with open(os.path.join(OUTDIR,'aggrdata','matrix.mtx'),'rb') as mtx_in:
		with gzip.open(os.path.join(OUTDIR,'aggrdata','matrix.mtx.gz'),'wb') as mtx_gz:
			shutil.copyfileobj(mtx_in, mtx_gz)
	os.remove(os.path.join(OUTDIR,'aggrdata','matrix.mtx')) 

	# write barcodes
	with gzip.open(os.path.join(OUTDIR,'aggrdata','barcodes.tsv.gz'), 'wb') as barcode_file:
		for barcode in BCS:
			barcode_file.write('{}\n'.format(barcode).encode())
	        
	# copy features file
	shutil.copyfile(os.path.join(config['INDIVDIR'], samples[0], config['SUBPATH_IN_INDIVDIRS'], 'features.tsv.gz'), os.path.join(OUTDIR,'aggrdata','features.tsv.gz'))

def aggrmatrix(config):
	OUTDIR=os.path.join(os.path.abspath(config['OUTDIR']),'')
	samples = next(os.walk(config['INDIVDIR']))[1]

	tocheck = []
	tocheck.append(os.path.join(OUTDIR,'aggrdata','barcodes.tsv.gz'))
	tocheck.append(os.path.join(OUTDIR,'aggrdata','features.tsv.gz'))
	tocheck.append(os.path.join(OUTDIR,'aggrdata','matrix.mtx.gz'))
	tocheck = [os.path.exists(f) for f in tocheck]
	allexist = all(v==True for v in tocheck)

	if allexist==False:
		compute_aggregation(config,samples)
	else:
		pass

	return samples