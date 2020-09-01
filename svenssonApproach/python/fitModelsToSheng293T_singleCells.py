import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

# Removes length column
datatypes = ["exonR", "exonU", "intronR", "intronU", "totalR", "totalU"]

dfs = [pd.read_csv('../Data/shengMatqseq/processedData/293T_clone1/%s.csv' % x, sep=',', index_col=0) for x in datatypes]

datasets = [anndata.AnnData(X=x.to_numpy().transpose(), var=pd.DataFrame(index=x.index), obs=pd.DataFrame(index=x.columns)) for x in dfs]
datasets = [x[:,np.sum(x.X,0)>0] for x in datasets]

for i,x in enumerate(datasets):	
	print("\n%s\n" % datatypes[i])
	x.uns['name'] = "Sheng et al. 2017"
	fitPoissonModel(x)	
	fitNBModels(x)
	x.obs['total_count'] = sum(x.X.transpose())
	x.write('../Data/output/sheng/293T_singleCell_%s.h5ad' % datatypes[i])

makePlotOfPoissonAndNB(datasets, datatypes, "../Figures/shengSingleCells293T.pdf")
