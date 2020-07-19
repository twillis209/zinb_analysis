import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *	
# Originally thought about adding the stats to existing files, but perhaps we should just use the output of my R code?
kolodAD=anndata.read('../Data/output/kolodERCC.h5ad')
ziDF=pd.read_csv('../Data/output/kolodERCCzi.csv', index_col=0)

kolodAD=kolodAD[:, ziDF.index]

# ml_mean differs for NB and ZINB models
kolodAD.var=pd.merge(kolodAD.var, ziDF['genewise_zi_zero_fraction'], left_index=True, right_index=True)
