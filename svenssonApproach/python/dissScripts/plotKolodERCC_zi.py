import sys
import os
import anndata

sys.path.append('..')
from svenssonCode import *

ad = anndata.read('../../Data/output/kolodZi/kolodERCC_zi.h5ad')
ad.uns['name'] = "Kolodziejczyk et al. 2015 (Smart-seq)"

makePlotOfNBWithZi([ad], ["ERCC spike-ins"], "../../Figures/kolodERCC_zi.pdf", figSize=(12, 8), nrows=1, ncols=1, title=None, limits=(1e-2, 1e6))
