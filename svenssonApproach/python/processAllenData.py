import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

allenDF=pd.read_csv('../Data/allenTophatCounts.csv', sep=',', index_col=0)
allenMetadataDF=pd.read_csv('../Data/allenMetadata.csv', sep=',', index_col=0)

allenAD=anndata.AnnData(X=allenDF.to_numpy().transpose(), var=pd.DataFrame(index=allenDF.index), obs=allenMetadataDF)

allenAD.write('../Data/output/allenData.h5ad')
