from svenssonCode import *
import os

filenames = [
"patel_CSC6.h5ad",
"patel_CSC8.h5ad",
"patel_MGH26-2.h5ad",
"patel_MGH26.h5ad",
"patel_MGH28.h5ad",
"patel_MGH29.h5ad",
"patel_MGH30.h5ad",
"patel_MGH30L.h5ad",
"patel_MGH31.h5ad"
]

datasets = [anndata.read(os.path.join('../Data/output', x)) for x in filenames]

for x in datasets:
	x.uns['name'] = 'Patel et al. 2014'

annotations = [
"CSC6",
"CSC8",
"MGH26 (2)",
"MGH26",
"MGH28",
"MGH29",
"MGH30",
"MGH30L",
"MGH31"
]

makePlotOfNBFits(datasets, annotations, outputPath="patelFits.pdf", nrows=3, ncols=3, figSize=(40,25))
