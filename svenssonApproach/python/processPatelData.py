import anndata
import re
import pandas as pd
from svenssonCode import *

patelDF=pd.read_csv('../Data/glioblastoma_raw_rnaseq_SCandbulk_counts_withannots.txt', sep="\t", index_col=2)

sraDF=pd.read_csv('../Data/glioblastomaRunTable.txt', sep=',', index_col=0)

# Metadata contains a source_name column, but the labels in the count matrix actually have more information wrt. to batch in the case of MGH26
patelAD=anndata.AnnData(X=patelDF.iloc[1:, 2:].to_numpy().transpose(), var=pd.DataFrame(data={"Symbol":patelDF.iloc[1:,1]}), obs=sraDF[["patient_id", "subtype", "Cell_Line", "AvgSpotLen"]])

patelAD.obs['source_name']=patelDF.iloc[0,2:]

# Filter from Risso
#patelAD=patelAD[:, np.sum(patelAD.X>50, axis=0)>=50]

substrings = ['MGH26_', 'MGH26-2_', 'MGH28_', 'MGH29_', 'MGH30_', 'MGH30L_', 'MGH31_', 'CSC6_', 'CSC8_']

patelIndices=[patelAD.obs.index[[i for i,x in enumerate(patelAD.obs.source_name) if re.search(y, x)]] for y in substrings]

datasets=[patelAD[x,:] for x in patelIndices]

for x in zip(datasets, [y[:-1] for y in substrings]):
	x[0].uns['shortName'] = x[1]
	x[0].uns['name'] = 'Patel et al. 2014 (%s)' % x[1]
#	fitPoissonModel(x[0])
#w	fitNBModels(x[0])
	x[0].write('../Data/output/patel_%s.h5ad' % x[1])
