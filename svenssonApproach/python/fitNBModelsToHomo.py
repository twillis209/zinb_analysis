import anndata
import pandas as pd
import numpy as np
import re 
from fitNBModels import *

kolodDF = pd.read_csv('../Data/kolodziejczyk_counttable_es.csv', sep=' ', index_col=0)
kolodAD=anndata.AnnData(X=kolodDF.to_numpy().transpose(), var=pd.DataFrame(index=kolodDF.index), obs=pd.DataFrame(index=kolodDF.columns))

prefixes=['ola_mES_lif', 'ola_mES_2i', 'ola_mES_a2i', 'ola_mES_lif_1', 'ola_mES_lif_2', 'ola_mES_lif_3', 'ola_mES_2i_2', 'ola_mES_2i_3', 'ola_mES_2i_4','ola_mES_2i_5', 'ola_mES_a2i_2', 'ola_mES_a2i_3']

# TODO: Need to add ERCC index to this, too
kolodIndices=[kolodAD.obs.index[[i for i,x in enumerate(kolodAD.obs.index) if re.search(y, x)]] for y in batchPrefixes]

#kolodLif=kolodAD[lifIndex, :]
#kolodi2=kolodAD[i2Index, :]
#koloda2i=kolodAD[a2iIndex, :]
#
#kolodLif.uns['name'] = 'Kolodziejczyk et al. 2015 (2)'
#kolodLif.uns['shortName'] = 'kolodLifHomo'
#kolodi2.uns['name'] = 'Kolodziejczyk et al. 2015 (3)'
#kolodi2.uns['shortName'] = 'kolod2iHomo'
#koloda2i.uns['name'] = 'Kolodziejczyk et al. 2015 (4)'
#koloda2i.uns['shortName'] = 'koloda2iHomo'

datasets=[kolodAD[x,:] for x in kolodIndices]

for x in zip(datasets, [re.sub('^ola_mES_', '', y) for y in batchPrefixes]):
	x[0].uns['shortName'] = x[1]
	x[0].uns['name'] = 'Kolodziejczyk et al. 2015 (%s)' % x[1]
	fitNBModels(x[0])
	anndata.write('../Data/output/kolod_%s.h5ad' % x[1])
