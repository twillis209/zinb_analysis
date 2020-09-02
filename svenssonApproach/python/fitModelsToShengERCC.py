import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

# Removes length column
shengDF = pd.read_csv('../Data/shengMatqseq/processedData/allERCC.csv', sep=',', index_col=0).iloc[:, 1:]

# Calculated using DESeq2's 'estimateSizeFactorsForMatrix' function, not using these at the moment
sizeFactors = [0.7648592, 1.1928399, 1.1955849, 1.3242736, 1.0770043, 0.5951686]

shengAD=anndata.AnnData(X=shengDF.to_numpy().transpose(), var=pd.DataFrame(index=shengDF.index), obs=pd.DataFrame(index=shengDF.columns))

shengAD.uns['name'] = "Sheng et al. 2017"

fitPoissonModel(shengAD)

fitNBModels(shengAD)

shengAD.obs['total_count'] = sum(shengAD.X.transpose())

shengAD.obs['size_factors'] = sizeFactors

shengAD.write('../Data/output/sheng/ercc.h5ad')
