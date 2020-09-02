import sys
import anndata

sys.path.append('..')
from svenssonCode import *

shengAD=anndata.read('../../Data/output/sheng/ercc.h5ad')

makePlotOfPoissonAndNB([shengAD], ["ERCC spike-ins"], "../../Figures/shengERCCPlots.pdf", pLimits=(10e-9,10e-3), nbLimits=(1, 10e3))
