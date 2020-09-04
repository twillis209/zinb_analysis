import anndata
import re
import pandas as pd
from svenssonCode import *


substrings = ['MGH26', 'MGH26-2', 'MGH28', 'MGH29', 'MGH30', 'MGH30L', 'MGH31', 'CSC6', 'CSC8']

datasets=[anndata.read('../Data/output/patel_%s.h5ad' % x) for x in substrings]

for i,x in enumerate(datasets):
	fitPoissonModel(x)
	fitNBModels(x)
	x.write('../Data/output/patel_%s.h5ad' % substrings[i])
