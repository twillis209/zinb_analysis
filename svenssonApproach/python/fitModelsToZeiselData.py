import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

#Categories (7, object): [astrocytes_ependymal, endothelial-mural, interneurons, microglia, oligodendrocytes,
#                         pyramidal CA1, pyramidal SS]

zeiselAD = anndata.read('../Data/output/zeiselData.h5ad')

vsmcAD = zeiselAD[zeiselAD.obs.loc[zeiselAD.obs.tissue=="Vsmc"].index]

epenAD = zeiselAD[zeiselAD.obs.loc[zeiselAD.obs.cell_type=="astrocytes_ependymal"].index]

for x in zip(['vsmc', 'epen'], [vsmcAD, epenAD]):
	fitNBModels(x[1])
	x[1].write('../Data/output/zeisel_%s.h5ad' % x[0])
