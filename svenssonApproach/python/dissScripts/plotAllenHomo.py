
import os
import sys

sys.path.append('..')

from svenssonCode import *

scnn1aAD = anndata.read('../../Data/output/allen_scnn1a.h5ad')

scnn1aAD.uns['name'] = 'Tasic et al. 2016'

makePlotOfNBFits([scnn1aAD], ["L4 Scnn1a cells"], outputPath="../../Figures/allenHomo.pdf", nrows=1, ncols=1, figSize=(10,8), limits=(10e-3, 10e3))
