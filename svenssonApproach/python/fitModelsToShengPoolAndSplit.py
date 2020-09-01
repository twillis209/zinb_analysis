import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

datatypes = ["exonR", "exonU", "intronR", "intronU", "totalR", "totalU"]

dfs = [pd.read_csv('../Data/shengMatqseq/processedData/setOne293T_sAndp/setOne293T_%s.csv' % x, sep=',', index_col=0) for x in datatypes]

datasets = [anndata.AnnData(X=x.to_numpy().transpose(), var=pd.DataFrame(index=x.index), obs=pd.DataFrame(index=x.columns)) for x in dfs]
datasets = [x[:,np.sum(x.X,0)>0] for x in datasets]

for i,x in enumerate(datasets):	
	print("\n%s\n" % datatypes[i])
	x.uns['name'] = "Sheng et al. 2017"
	fitPoissonModel(x)	
	fitNBModels(x)
	x.obs['total_count'] = sum(x.X.transpose())
	x.write('../Data/output/sheng/setOne293T_sAndp_%s.h5ad' % datatypes[i])

makePlotOfPoissonAndNB(datasets, datatypes, "../Figures/shengSetOne293T.pdf")

dfs = [pd.read_csv('../Data/shengMatqseq/processedData/setTwo293T_sAndp/setTwo293T_%s.csv' % x, sep=',', index_col=0) for x in datatypes]

datasets = [anndata.AnnData(X=x.to_numpy().transpose(), var=pd.DataFrame(index=x.index), obs=pd.DataFrame(index=x.columns)) for x in dfs]
datasets = [x[:,np.sum(x.X,0)>0] for x in datasets]

for i,x in enumerate(datasets):	
	print("\n%s\n" % datatypes[i])
	x.uns['name'] = "Sheng et al. 2017"
	fitPoissonModel(x)	
	fitNBModels(x)
	x.obs['total_count'] = sum(x.X.transpose())
	x.write('../Data/output/sheng/setTwo293T_sAndp_%s.h5ad' % datatypes[i])

makePlotOfPoissonAndNB(datasets, datatypes, "../Figures/shengSetTwo293T.pdf")
