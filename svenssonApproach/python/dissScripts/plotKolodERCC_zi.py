import sys
import os
import anndata

sys.path.append('..')
from svenssonCode import *

ad = anndata.read('../../Data/output/kolodZi/kolodERCC_zi.h5ad')

makePlotOfNBWithZi([ad], ["ERCC spike-ins"], "../../Figures/kolodERCC_zi.pdf", figSize=(15, 8), nrows=1, ncols=1, title=None, limits=(10e-2, 10e5))
