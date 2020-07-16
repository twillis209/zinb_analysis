'''
Just experimentation at the moment.
'''
import statsmodels.api as sm
import numpy as np
import pandas as pd
import anndata

fletcherAD=anndata.read('svenssonApproach/Data/output/fletcherERCC.h5ad')
testDF=pd.DataFrame(fletcherAD.X[0,:], columns=['Counts'])
testModeD=sm.ZeroInflatedNegativeBinomialP.from_formula(formula='Counts ~ 1', data=testDF)
