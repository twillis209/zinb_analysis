import os
import sys

sys.path.append('..')
from svenssonCode import *

filenames = [
	"kolod_2i_2_zi.h5ad"
#	"kolod_2i_zi.h5ad",
#	"kolodERCC_zi.h5ad",
#	"kolodNoERCC_zi.h5ad"	
		]

datasets = [anndata.read(os.path.join('../../Data/output/kolodZi', x)) for x in filenames]

	#
annotations = [
	'\'2i\' medium, batch 1'
]

datasets[0].uns['name'] =  'Kolodziejczyk et al. 2015'

makePlotOfNBWithZi(datasets, annotations, outputPath="../../Figures/kolodHomoZi.pdf", nrows=1, ncols=1, figSize=(15,8), limits=(10e-3, 10e3))
