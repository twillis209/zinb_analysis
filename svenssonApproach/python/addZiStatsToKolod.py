import anndata
import pandas as pd
import numpy as np
import re 
import os
from svenssonCode import *	
# ../Data/output/*

h5adFilenames = [
		"kolod_2i_2.h5ad",
		"kolod_2i_3.h5ad",
		"kolod_2i_4.h5ad",
		"kolod_2i_5.h5ad",
		"kolod_2i.h5ad",
		"kolod_a2i_2.h5ad",
		"kolod_a2i_3.h5ad",
		"kolod_a2i.h5ad",
		"kolod_lif_1.h5ad",
		"kolod_lif_2.h5ad",
		"kolod_lif_3.h5ad",
		"kolod_lif.h5ad",
		"kolodERCC.h5ad"
	]

csvFilenames = [
		"kolodZi_ola_mES_2i_2.csv",
		"kolodZi_ola_mES_2i_3.csv",
		"kolodZi_ola_mES_2i_4.csv",
		"kolodZi_ola_mES_2i_5.csv",
		"kolodZi_ola_mES_2i.csv",
		"kolodZi_ola_mES_a2i_2.csv",
		"kolodZi_ola_mES_a2i_3.csv",
		"kolodZi_ola_mES_a2i.csv",
		"kolodZi_ola_mES_lif_1.csv",
		"kolodZi_ola_mES_lif_2.csv",
		"kolodZi_ola_mES_lif_3.csv",
		"kolodZi_ola_mES_lif.csv",
		"kolodERCCzi.csv"
		]

# ../Data/output/kolodZi*

# ml_mean differs for NB and ZINB models
for pair in zip(h5adFilenames, csvFilenames):
	ad = anndata.read(os.path.join("../Data/output", pair[0]))
	ziDF = pd.read_csv(os.path.join("../Data/output/kolodZi", pair[1]), index_col=0)
	ad = ad[:, ziDF.index]
	ad.var=pd.merge(ad.var, ziDF['genewise_zi_zero_fraction'], left_index=True, right_index=True)
	ad.uns['shortName'] = pair[0].replace('.h5ad', '')
	ad.write(os.path.join("../Data/output/kolodZi", ad.uns['shortName']+'_zi.h5ad'))
