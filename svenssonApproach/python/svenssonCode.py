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
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

"""
Code for fitting negative binomial models and plotting the results. Adapted from Svensson 2020.
"""

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

def makePlot(datasets, annotationDict, locationDict, outputPath, figSize=(40, 25), nrows=4, ncols=4):
	mins = []
	maxs = []

	for adata in datasets:
		difference1 = adata.var['empirical_zero_fraction'] - adata.var['global_zero_fraction']
		difference2 = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
		difference3 = adata.var['empirical_zero_fraction'] - adata.var['genewise_zi_zero_fraction']
		mins.append(min(difference1.min(), difference2.min(), difference3.min()))
		maxs.append(max(difference1.max(), difference2.max(), difference3.max()))

	min(mins), max(maxs)
	
	# TODO: still editing size and grid dimensions by hand
	fig = plt.figure(figsize=figSize)

	outer_grid = fig.add_gridspec(nrows=nrows, ncols=ncols, hspace=0.3, wspace=0.2)

	for adata in datasets:
	    
		i = locationDict[adata.uns['name']]
		grid_box = outer_grid[i]

		inner_grid = grid_box.subgridspec(2, 3)
		# 0 1 2
		# 3 4 5 	
		# -- Global -- 

		ax = fig.add_subplot(inner_grid[0])

		ax.set_title('Common dispersion')

		ax.set_xscale('log')
		ax.set_xlim(left=2e-4, right=1e4)

		ax.scatter(adata.var['empirical_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', label='Observed', rasterized=True);
		ax.scatter(adata.var['empirical_mean'],
		       adata.var['global_zero_fraction'],
		       ec='w', c='grey', label='Expected', rasterized=True);

		ax.set_ylabel('Fraction zeros')

		ax.legend(title='Genes', loc='lower left', scatterpoints=3, fontsize=8)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Global dispersion differences

		ax = fig.add_subplot(inner_grid[3])

		ax.set_xscale('log')
		ax.set_xlim(left=2e-4, right=1e4)
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['global_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', label='Genes', rasterized=True)

		ax.set_ylabel('Difference \n(Observed - Expected)')
		ax.set_xlabel('Mean')

		ax.legend(loc='lower left', scatterpoints=3)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)


		# -- Genewise --

		ax = fig.add_subplot(inner_grid[1])

		ax.set_title('Gene-wise dispersion')

		ax.set_xscale('log')
		ax.set_xlim(left=2e-4, right=1e4)

		ax.scatter(adata.var['empirical_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', rasterized=True);
		ax.scatter(adata.var['empirical_mean'],
		       adata.var['genewise_zero_fraction'],
		       ec='w', c='grey', rasterized=True);

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Genewise dispersion differences

		ax = fig.add_subplot(inner_grid[4])

		ax.set_xscale('log')
		ax.set_xlim(left=2e-4, right=1e4)
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', rasterized=True)

		ax.set_xlabel('Mean')

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Zero inflation

		ax = fig.add_subplot(inner_grid[2])

		ax.set_title('ZI model')

		ax.set_xscale('log')
		ax.set_xlim(left=2e-4, right=1e4)

		ax.scatter(adata.var['empirical_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', rasterized=True);
		ax.scatter(adata.var['empirical_mean'],
		       adata.var['genewise_zi_zero_fraction'],
		       ec='w', c='grey', rasterized=True);

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Zero inflation differences

		ax = fig.add_subplot(inner_grid[5])

		ax.set_xscale('log')
		ax.set_xlim(left=2e-4, right=1e4)
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['genewise_zi_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', rasterized=True)

		ax.set_xlabel('Mean')

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		# -- Annotation --

		bbox = grid_box.get_position(fig)
		x = (bbox.x0 + bbox.x1) / 2
		y = bbox.y1 + 0.015

		fig.text(x, y, adata.uns['name'] + '\n' + annotationDict[adata.uns['name']], ha='center', fontsize=12)
	    
	fig.savefig(outputPath, dpi=500, bbox_inches='tight')
