import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

kolodAD=anndata.read('../Data/output/kolodNoERCC.h5ad')


kolodAD.uns['shortName'] = 'kolodNoERCC'
kolodAD.uns['name'] = 'Kolodziejczyk et al. 2015 (all bio. features)'
fitNBModels(kolodAD)
kolodAD.write('../Data/output/kolodNoERCC.h5ad')
