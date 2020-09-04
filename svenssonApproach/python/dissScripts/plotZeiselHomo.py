
import os
import sys

sys.path.append('..')

from svenssonCode import *

epenAD = anndata.read('../../Data/output/zeisel_epen.h5ad')
vsmcAD = anndata.read('../../Data/output/zeisel_vsmc.h5ad')

epenAD.uns['name'] = 'Zeisel et al. 2015'
vsmcAD.uns['name'] = 'Zeisel et al. 2015'

#makePlotOfNBFits([epenAD, vsmcAD], ["VSMC", "Ependyma"], outputPath="../../Figures/zeiselHomo.pdf", nrows=2, ncols=1, figSize=(16,16), limits=(10e-3, 10e2))


makePlotOfNBFits([vsmcAD], ["Vascular smooth muscle cells"], outputPath="../../Figures/zeiselHomo.pdf", nrows=1, ncols=1, figSize=(10,8), limits=(10e-3, 10e2))
