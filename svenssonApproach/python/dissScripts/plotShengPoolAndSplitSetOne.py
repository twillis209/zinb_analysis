import sys
import anndata
import pandas as pd
import numpy as np
import re 


sys.path.append('..')
from svenssonCode import *

datasets = [anndata.read('../../Data/output/sheng/setOne293T_sAndp_%s.h5ad' % x) for x in ['exonR', 'exonU']]

annotations = ["Read counts (exonic)", "UMI counts (exonic)"]

makePlotOfNBFits(datasets, annotations, "../../Figures/shengSetOne293T.pdf", title="Pool-and-split HEK293T samples, Sheng et al. 2017", figSize=(30,25))

datasets = [anndata.read('../../Data/output/sheng/noFiltersetOne293T_sAndp_%s.h5ad' % x) for x in ['exonR', 'exonU']]

makePlotOfNBFits(datasets, annotations, "../../Figures/noFilterShengSetOne293T.pdf", title="Pool-and-split HEK293T samples, Sheng et al. 2017", figSize=(30,25))
