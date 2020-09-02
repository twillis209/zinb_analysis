import sys
import anndata

sys.path.append('..')
from svenssonCode import *

datasets = [anndata.read('../../Data/output/sheng/293T_singleCell_%s.h5ad' % x) for x in ['exonR', 'exonU']]

annotations = ["Read counts (exonic)", "UMI counts (exonic)"]

makePlotOfNBFits(datasets, annotations, "../../Figures/sheng293T_singleCell.pdf", title="HEK293T single cells, Sheng et al. 2017", figSize=(30,25))

datasets = [anndata.read('../../Data/output/sheng/noFilter293T_singleCell_%s.h5ad' % x) for x in ['exonR', 'exonU']]

makePlotOfNBFits(datasets, annotations, "../../Figures/noFilterSheng293T_singleCell.pdf", title="HEK293T single cells, Sheng et al. 2017", figSize=(30,25))
