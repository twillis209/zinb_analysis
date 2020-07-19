from glob import glob
from svenssonCode import *
import anndata
import pandas as pd
import matplotlib.pyplot as plt

datasets = [
	anndata.read('../Data/output/allenERCC.h5ad'),
	anndata.read('../Data/output/fletcherERCC.h5ad'),
	anndata.read('../Data/output/kolodERCC.h5ad'),
	anndata.read('../Data/output/zhengERCC.h5ad'),
	anndata.read('../Data/output/zeiselERCC.h5ad'),
	anndata.read('../Data/output/kolod_2i_2.h5ad'),
	anndata.read('../Data/output/kolod_2i_3.h5ad'),
	anndata.read('../Data/output/kolod_2i_4.h5ad'),
	anndata.read('../Data/output/kolod_2i_5.h5ad'),
	anndata.read('../Data/output/kolod_2i.h5ad'),
	anndata.read('../Data/output/kolod_a2i_2.h5ad'),
	anndata.read('../Data/output/kolod_a2i_3.h5ad'),
	anndata.read('../Data/output/kolod_a2i.h5ad'),
	anndata.read('../Data/output/kolod_lif_1.h5ad'),
	anndata.read('../Data/output/kolod_lif_2.h5ad'),
	anndata.read('../Data/output/kolod_lif_3.h5ad'),
	anndata.read('../Data/output/kolod_lif.h5ad')
]

names = [ 
	'Tasic et al. 2016',
	'Fletcher et al. 2017',
	'Kolodziejczyk et al. 2015 (ERCC)',
	'Zheng et al. 2017',
	'Zeisel et al. 2015',
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
	'Kolodziejczyk et al. 2015 (lif)'
]

for i,x in enumerate(datasets):
	x.uns['name'] = names[i]

annotation = {
	'Tasic et al. 2016' : 'ERCC spike-ins',
	'Fletcher et al. 2017' :'ERCC spike-ins', 
	'Kolodziejczyk et al. 2015 (ERCC)' :'ERCC spike-ins', 
	'Zheng et al. 2017' : 'ERCC spike-ins',
	'Zeisel et al. 2015' : 'ERCC spike-ins',
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

location = {
	'Tasic et al. 2016':0,
	'Fletcher et al. 2017':1,
	'Kolodziejczyk et al. 2015 (ERCC)':2,
	'Zheng et al. 2017':3,
	'Zeisel et al. 2015':4,
	'Kolodziejczyk et al. 2015 (2i_2)':5,
	'Kolodziejczyk et al. 2015 (2i_3)':6,
	'Kolodziejczyk et al. 2015 (2i_4)':7,
	'Kolodziejczyk et al. 2015 (2i_5)':8,
	'Kolodziejczyk et al. 2015 (2i)':9,
	'Kolodziejczyk et al. 2015 (a2i_2)':10,
	'Kolodziejczyk et al. 2015 (a2i_3)':11,
	'Kolodziejczyk et al. 2015 (a2i)':12,
	'Kolodziejczyk et al. 2015 (lif_1)':13,
	'Kolodziejczyk et al. 2015 (lif_2)':14,
	'Kolodziejczyk et al. 2015 (lif_3)':15,
	'Kolodziejczyk et al. 2015 (lif)':16
}

makePlot(datasets, annotation, location, outputPath='test.pdf')
