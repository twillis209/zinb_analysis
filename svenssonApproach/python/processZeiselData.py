import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

zeiselDF=pd.read_csv('../Data/zeiselCounts.csv', sep=',', index_col=0)
zeiselCellTypesDF=pd.read_csv('../Data/zeiselCellType.csv', sep=',', index_col=0)

zeiselAD=anndata.AnnData(X=zeiselDF.to_numpy().transpose(), var=pd.DataFrame(index=zeiselDF.index), obs=pd.DataFrame(index=zeiselDF.columns))

zeiselCellTypesDF.index = zeiselDF.columns

zeiselAD.obs = zeiselCellTypesDF

zeiselAD.write('../Data/output/zeiselData.h5ad')
