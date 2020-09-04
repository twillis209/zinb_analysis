import sys
import os
import anndata

sys.path.append('..')
from svenssonCode import *

inputRoot='../../Data/output'

filenames = [
	"allenERCC.h5ad",
	"fletcherERCC.h5ad",
	"kolodERCC.h5ad",
	"zeiselERCC.h5ad",
	"zhengERCC.h5ad"
]

datasets = [anndata.read(os.path.join(inputRoot, x)) for x in filenames]

allenAD=datasets[0]
fletcherAD=datasets[1]
kolodAD=datasets[2]
zeiselAD=datasets[3]
zhengAD=datasets[4]

#makePlotOfPoissonFits(datasets, ["ERCC spike-ins"]*5, "../../Figures/rissoPoissonERCCPlots.pdf", nrows=2, ncols=3, figSize=(20	, 25), pLimits=(10e-10, 1))
#
#makePlotOfNBFits(datasets, ["ERCC spike-ins"]*5, "../../Figures/rissoNBERCCPlots.pdf", nrows=3, ncols=2, figSize=(20, 25))

makePlotOfNBFits([allenAD, fletcherAD], ["ERCC spike-ins"]*2, "../../Figures/rissoNBERCCPlots_allenFletcherKolod.pdf", nrows=1, ncols=2, figSize=(20, 12), limits=(10e-4, 10e4))

#makePlotOfPoissonFits([zeiselAD, zhengAD], ["ERCC spike-ins"]*2, "../../Figures/rissoPoissonERCCPlots_zeiselZheng.pdf", nrows=1, ncols=2, figSize=(10, 15), pLimits=(10e-10, 1))
