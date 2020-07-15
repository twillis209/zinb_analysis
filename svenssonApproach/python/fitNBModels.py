import anndata
import pandas as pd
import numpy as np
import re 
from scipy import optimize
from scipy.special import gammaln
from scipy.special import psi
from scipy.special import factorial
from scipy.optimize import fmin_l_bfgs_b as optim
from tqdm import tqdm

# X is a numpy array representing the data
# initial params is a numpy array representing the initial values of
# size and prob parameters
def fit_nbinom(X, initial_params=None):
    ''' This code is adapted from https://github.com/gokceneraslan/fit_nbinom
    '''
    infinitesimal = np.finfo(np.float).eps

    def log_likelihood(params, *args):
        r, p = params
        X = args[0]
        N = X.size

        # MLE estimate based on the formula on Wikipedia:
        # http://en.wikipedia.org/wiki/Negative_binomial_distribution#Maximum_likelihood_estimation
        result = np.sum(gammaln(X + r)) - np.sum(np.log(factorial(X)))- N * (gammaln(r))+ N * r * np.log(p)+ np.sum(X * np.log(1 - (p if p < 1 else 1 - infinitesimal)))

        return -result

    if initial_params is None:
        # reasonable initial values (from fitdistr function in R)
        m = np.mean(X)
        v = np.var(X)
        size = (m ** 2) / (v-m) if v > m else 10

        # convert mu/size parameterization to prob/size
        p0 = size / ((size + m) if size + m != 0 else 1)
        r0 = size
        initial_params = np.array([r0, p0])

    bounds = [(infinitesimal, None), (infinitesimal, 1)]
    optimres = optim(log_likelihood,
                     x0=initial_params,
                     args=(X,),
                     approx_grad=1,
                     bounds=bounds)

    params = optimres[0]
    return {'size': params[0], 'prob': params[1]}


# First determine the empirical statistics for each gene.
# 
# Then fit NB parameters to the each gene by ML.
def fit_per_gene_stats(adata):
    adata.var['empirical_mean'] = np.nan
    adata.var['empirical_variance'] = np.nan
    adata.var['empirical_zero_fraction'] = np.nan
    adata.var['ml_mean'] = np.nan
    adata.var['genewise_dispersion'] = np.nan

    for g in tqdm(adata.var.index):
        y = adata[:, g].X
        adata.var.loc[g, 'empirical_mean'] = y.mean(0)
        adata.var.loc[g, 'empirical_variance'] = y.var(0)
        adata.var.loc[g, 'empirical_zero_fraction'] = 1 - (y > 0).sum(0) / y.shape[0]
        result = fit_nbinom(y)
        mu = ((1. - result['prob']) * result['size']) / (result['prob'])
        adata.var.loc[g, 'ml_mean'] = mu
        adata.var.loc[g, 'genewise_dispersion'] = 1. / result['size']


# Now fit global parameters to each data set

def var_fun(mu, phi):
    return mu + phi * mu ** 2

def prob_zero_fun(mu, phi):
    if phi == .0:
        return np.exp(-mu)
    
    phi_1 = 1. / phi
    return (phi_1 / (mu + phi_1)) ** phi_1

def prob_zero_fun_vec(mu, phi):
    phi_1 = 1. / phi
    return (phi_1 / (mu + phi_1)) ** phi_1

def fitNBModels(adata):
	fit_per_gene_stats(adata)
	phi_hat, _ = optimize.curve_fit(var_fun, adata.var['empirical_mean'], adata.var['empirical_variance'])
	adata.uns['global_dispersion'] = phi_hat
	adata.var['global_zero_fraction'] = prob_zero_fun_vec(adata.var['empirical_mean'],
							  adata.uns['global_dispersion'])
	adata.var['genewise_zero_fraction'] = prob_zero_fun_vec(adata.var['empirical_mean'],
							    adata.var['genewise_dispersion'])
	return
