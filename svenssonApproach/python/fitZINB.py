'''
Just experimentation at the moment.
'''
import statsmodels.api as sm
import numpy as np
import pandas as pd
import anndata

fletcherAD=anndata.read('../Data/output/fletcherERCC.h5ad')
testDF=pd.DataFrame(fletcherAD.X[2,:], columns=['Counts'])
endog=testDF.to_numpy()
exog=np.array([1]*92)
testModel=sm.ZeroInflatedNegativeBinomialP(endog=endog, exog=exog)
testModel.fit()
