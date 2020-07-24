from svenssonCode import *
import os


filenames = [
	"kolod_2i_2_zi.h5ad",
	"kolod_2i_3_zi.h5ad",
	"kolod_2i_4_zi.h5ad",
	"kolod_2i_5_zi.h5ad",
	"kolod_2i_zi.h5ad",
	"kolod_a2i_2_zi.h5ad",
	"kolod_a2i_3_zi.h5ad",
	"kolod_a2i_zi.h5ad",
	"kolod_lif_1_zi.h5ad",
	"kolod_lif_2_zi.h5ad",
	"kolod_lif_3_zi.h5ad",
	"kolod_lif_zi.h5ad",
	"kolodERCC_zi.h5ad",
	"kolodNoERCC_zi.h5ad"
		]

filenames = [
	"kolod_2i_2_zi.h5ad",
	"kolod_2i_zi.h5ad",
	"kolodERCC_zi.h5ad",
	"kolodNoERCC_zi.h5ad"	
		]

datasets = [anndata.read(os.path.join('../Data/output/kolodZi', x)) for x in filenames]

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
