import sys
import anndata

sys.path.append('..')
from svenssonCode import *

exonRad = anndata.read('../../Data/output/sheng/293T_singleCell_exonR.h5ad')
exonUad = anndata.read('../../Data/output/sheng/293T_singleCell_exonU.h5ad')

annotations = ["Read counts (exonic)", "UMI counts (exonic)"]

exonRad.uns['name'] = "Sheng et al. 2017 (MATQ-seq)"
exonUad.uns['name'] = "Sheng et al. 2017 (MATQ-seq)"

makePlotOfNBFits([exonRad], ["HEK293T read counts (exonic)"], outputPath="../../Figures/sheng293T_singleCell_exonR.pdf", figSize=(8,8), limits=(1e-3, 1e4), nrows=1, ncols=1)
makePlotOfNBFits([exonUad], ["HEK293T UMI counts (exonic)"], outputPath="../../Figures/sheng293T_singleCell_exonU.pdf",figSize=(8,8), limits=(1e-3, 1e4), nrows=1, ncols=1)
