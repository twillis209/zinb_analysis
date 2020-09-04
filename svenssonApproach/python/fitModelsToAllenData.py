import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

allenAD = anndata.read('../Data/output/allenData.h5ad')

scnn1aAD = allenAD[allenAD.obs.loc[(allenAD.obs['Primary.Type'] == 'L4 Scnn1a') & (allenAD.obs['Core.Type'] == 'Core') & (allenAD.obs.dissection_s == 'L4')].index]

scnn1aAD = scnn1aAD[:,[i for i,x in enumerate(scnn1aAD.var.index) if not re.search('^ERCC', x)]]

fitNBModels(scnn1aAD)

scnn1aAD.write('../Data/output/allen_scnn1a.h5ad')
