from svenssonCode import *
import os
	x[0].write('../Data/output/patel_%s.h5ad' % x[1])

filenames = [
	'patel_MGH26.h5ad',
	'patel_MGH26-2.h5ad',
	'patel_MGH28.h5ad',
	'patel_MGH29.h5ad',
	'patel_MGH30.h5ad',
	'patel_MGH30L.h5ad',
	'patel_MGH31.h5ad',
	'patel_CSC6.h5ad',
	'patel_CSC8.h5ad'
		]

datasets = [anndata.read(os.path.join('../Data/output/patel', x)) for x in filenames]

names = [ 
	'Kolodziejczyk et al. 2015 (2i_2)',
	'Kolodziejczyk et al. 2015 (2i_3)',
	'Kolodziejczyk et al. 2015 (2i_4)',
	'Kolodziejczyk et al. 2015 (2i_5)',
	'Kolodziejczyk et al. 2015 (2i)',
	'Kolodziejczyk et al. 2015 (a2i_2)',
	'Kolodziejczyk et al. 2015 (a2i_3)',
	'Kolodziejczyk et al. 2015 (a2i)',
	'Kolodziejczyk et al. 2015 (lif_1)',
	'Kolodziejczyk et al. 2015 (lif_2)',
	'Kolodziejczyk et al. 2015 (lif_3)',
	'Kolodziejczyk et al. 2015 (lif)',
	'Kolodziejczyk et al. 2015 (ERCC)',
	'Kolodziejczyk et al. 2015 (No ERCC)'
	]

names = [
	'Kolodziejczyk et al. 2015 (2i_2)',
	'Kolodziejczyk et al. 2015 (2i)',
	'Kolodziejczyk et al. 2015 (ERCC)',
	'Kolodziejczyk et al. 2015 (No ERCC)'
	]

for i,x in enumerate(datasets):
	x.uns['name'] = names[i]

annotation = {
	'Kolodziejczyk et al. 2015 (No ERCC)' :'All biological features', 
	'Kolodziejczyk et al. 2015 (ERCC)' :'ERCC spike-ins', 
	'Kolodziejczyk et al. 2015 (2i_2)' : '2i batch 1',
	'Kolodziejczyk et al. 2015 (2i_3)' : '2i batch 2',
	'Kolodziejczyk et al. 2015 (2i_4)' : '2i batch 3',
	'Kolodziejczyk et al. 2015 (2i_5)' : '2i batch 4', 
	'Kolodziejczyk et al. 2015 (2i)' : '2i medium',
	'Kolodziejczyk et al. 2015 (a2i_2)' : 'a2i batch 1',
	'Kolodziejczyk et al. 2015 (a2i_3)' : 'a2i batch 2',
	'Kolodziejczyk et al. 2015 (a2i)' : 'a2i medium',
	'Kolodziejczyk et al. 2015 (lif_1)' : 'lif batch 1',
	'Kolodziejczyk et al. 2015 (lif_2)' : 'lif batch 2',
	'Kolodziejczyk et al. 2015 (lif_3)' : 'lif batch 3',
	'Kolodziejczyk et al. 2015 (lif)' : 'lif medium'
}
annotation = {
	'Kolodziejczyk et al. 2015 (No ERCC)' :'All biological features', 
	'Kolodziejczyk et al. 2015 (ERCC)' :'ERCC spike-ins', 
	'Kolodziejczyk et al. 2015 (2i_2)' : '2i batch 1',
	'Kolodziejczyk et al. 2015 (2i)' : '2i medium',
}

location = dict(zip(names, range(len(names))))

makePlot(datasets, annotation, location, nrows=2, ncols=2, figSize=(20,13), outputPath='test.pdf')
