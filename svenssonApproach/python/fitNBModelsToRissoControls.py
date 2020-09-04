import anndata
import pandas as pd
import numpy as np
import re 
from svenssonCode import *

# Initialising AnnData objects using the guidance here: http://falexwolf.de/blog/171223_AnnData_indexing_views_HDF5-backing/
allenDF = pd.read_csv('../Data/allenTophatCounts.csv', sep=',', index_col=0)
allenAD=anndata.AnnData(X=allenDF.to_numpy().transpose(), var=pd.DataFrame(index=allenDF.index), obs=pd.DataFrame(index=allenDF.columns))
erccIndex=allenAD.var.index[[i for i,x in enumerate(allenAD.var.index) if re.search('ERCC', x)]]
allenAD=allenAD[:, erccIndex]

fletcherDF = pd.read_csv('../Data/GSE95601_oeHBCdiff_Cufflinks_eSet_counts_table.txt', sep='\t', index_col=0)
fletcherAD=anndata.AnnData(X=fletcherDF.to_numpy().transpose(), var=pd.DataFrame(index=fletcherDF.index), obs=pd.DataFrame(index=fletcherDF.columns))
erccIndex=fletcherAD.var.index[[i for i,x in enumerate(fletcherAD.var.index) if re.search('ERCC', x)]]
fletcherAD=fletcherAD[:, erccIndex]

kolodDF = pd.read_csv('../Data/kolodziejczyk_counttable_es.csv', sep=' ', index_col=0)
kolodAD=anndata.AnnData(X=kolodDF.to_numpy().transpose(), var=pd.DataFrame(index=kolodDF.index), obs=pd.DataFrame(index=kolodDF.columns))
erccIndex=kolodAD.var.index[[i for i,x in enumerate(kolodAD.var.index) if re.search('ERCC', x)]]
kolodAD=kolodAD[:, erccIndex]

# This only contains ERCC spike-in data
zeiselDF = pd.read_csv('../Data/zeiselSpikeins.csv', sep=',', index_col=0)
zeiselAD=anndata.AnnData(X=zeiselDF.to_numpy().transpose(), var=pd.DataFrame(index=zeiselDF.index), obs=pd.DataFrame(index=zeiselDF.columns))

zhengAD=anndata.read('../Data/zheng_gemcode_control.h5ad')

allenAD.uns['name'] = 'Tasic et al. 2016'
allenAD.uns['shortName'] = 'allenERCC'
fletcherAD.uns['name'] = 'Fletcher et al. 2017'
fletcherAD.uns['shortName'] = 'fletcherERCC'
kolodAD.uns['name'] = 'Kolodziejczyk et al. 2015 (1)'
kolodAD.uns['shortName'] = 'kolodERCC'
zhengAD.uns['name'] = 'Zheng et al. 2017'
zhengAD.uns['shortName'] = 'zhengERCC'
zeiselAD.uns['name'] = 'Zeisel et al. 2015'
zeiselAD.uns['shortName'] = 'zeiselERCC'

datasets=[allenAD, fletcherAD, kolodAD, zhengAD, zeiselAD]

for adata in datasets:
    fit_per_gene_stats(adata)

for adata in datasets:
    phi_hat, _ = optimize.curve_fit(var_fun, adata.var['empirical_mean'], adata.var['empirical_variance'])
    adata.uns['global_dispersion'] = phi_hat

## Calculate predicted values to compare the empirical values with
#
for adata in datasets:
    adata.var['global_zero_fraction'] = prob_zero_fun_vec(adata.var['empirical_mean'],
                                                          adata.uns['global_dispersion'])
    adata.var['genewise_zero_fraction'] = prob_zero_fun_vec(adata.var['empirical_mean'],
                                                            adata.var['genewise_dispersion'])
#
#
## Finally save the results for making plots
for adata in datasets:
    adata.write('../Data/output/' + adata.uns['shortName'] + '.h5ad')
