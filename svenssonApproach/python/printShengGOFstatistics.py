import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

datatypes = ["exonR", "exonU", "intronR", "intronU", "totalR", "totalU"]

for x in datatypes:
	print(x)
	for filePath in ['../Data/output/sheng/293T_singleCell_%s.h5ad', '../Data/output/sheng/setOne293T_sAndp_%s.h5ad', '../Data/output/sheng/setTwo293T_sAndp_%s.h5ad']:
		print(filePath % x)
		ad = anndata.read(filePath % x)
		printGOFstatistics(ad)	
