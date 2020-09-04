import os
import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

inputRoot='../Data/output'

filenames = [
        "allenERCC.h5ad",
        "fletcherERCC.h5ad",
        "kolodERCC.h5ad",
        "zeiselERCC.h5ad",
        "zhengERCC.h5ad"
]

datasets = [anndata.read(os.path.join(inputRoot, x)) for x in filenames]

for i,x in enumerate(datasets):
	fitPoissonModel(x)
	x.write(os.path.join(inputRoot, filenames[i]))
