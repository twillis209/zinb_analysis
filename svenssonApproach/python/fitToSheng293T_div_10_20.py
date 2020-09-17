import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

daf = pd.read_csv('../Data/shengMatqseq/processedData/293T_clone1/exonR.csv', sep=',', index_col=0)

ad = anndata.AnnData(X=daf.to_numpy().transpose(), var=pd.DataFrame(index=daf.index), obs=pd.DataFrame(index=daf.columns))

ad = ad[:,np.sum(ad.X,0)>0]

#for i,x in enumerate(datasets):	
#	print("\n%s\n" % datatypes[i])
#	x.uns['name'] = "Sheng et al. 2017"
#	fitPoissonModel(x)	
#	fitNBModels(x)
#	x.obs['total_count'] = sum(x.X.transpose())
#	x.write('../Data/output/sheng/293T_singleCell_%s.h5ad' % datatypes[i])
