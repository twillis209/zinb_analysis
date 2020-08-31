import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

# Removes length column
datatypes = ["exonR", "exonU", "intronR", "intronU", "totalR", "totalU"]

dfs = [pd.read_csv('../Data/shengMatqseq/processedData/setOne293T_%s.csv' % x, sep=',', index_col=0) for x in datatypes]

datasets = [anndata.AnnData(X=x.to_numpy().transpose(), var=pd.DataFrame(index=x.index), obs=pd.DataFrame(index=x.columns)) for x in dfs]

#fitPoissonModel(datasets[0])
fitNBModels(datasets[0])
#for i,x in enumerate(datasets):	
#	x.uns['name'] = "Sheng et al. 2017"
#	fitPoissonModel(x)	
#	fitNBModels(x)
#	x.obs['total_count'] = sum(x.X.transpose())
#
#makePlotOfPoissonAndNB(datasets, datatypes, "../Figures/shengSetOne293T.pdf")
