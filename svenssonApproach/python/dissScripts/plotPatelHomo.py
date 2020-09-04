
import os
import sys

sys.path.append('..')

from svenssonCode import *

ad = anndata.read('../../Data/output/patel_MGH26.h5ad')

ad.uns['name'] = 'Patel et al. 2014'

makePlotOfNBFits([ad], ["Glioblastoma cells"], outputPath="../../Figures/patelHomo.pdf", nrows=1, ncols=1, figSize=(8,8))
