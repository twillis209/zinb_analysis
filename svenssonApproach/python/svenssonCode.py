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

def fitPoissonModel(adata):
	def prob_zero_fun(mu, counts):
    		return np.exp(-(mu[:,np.newaxis]).dot(counts.values[np.newaxis])).sum(1) / len(counts)

	mean_of_scaled = np.array(adata.X.sum(0) / adata.X.sum()).reshape(-1)
	adata.var['scaled_count_mean'] = mean_of_scaled
	adata.obs['total_counts'] = np.array(adata.X.sum(1)).reshape(-1)
	adata.var['poisson_zero_fraction'] = prob_zero_fun(adata.var['scaled_count_mean'], adata.obs['total_counts'])
	return 

def fitNBModels(adata):
	fit_per_gene_stats(adata)
	phi_hat, _ = optimize.curve_fit(var_fun, adata.var['empirical_mean'], adata.var['empirical_variance'])
	adata.uns['global_dispersion'] = phi_hat
	adata.var['global_zero_fraction'] = prob_zero_fun_vec(adata.var['empirical_mean'],
							  adata.uns['global_dispersion'])
	adata.var['genewise_zero_fraction'] = prob_zero_fun_vec(adata.var['empirical_mean'],
							    adata.var['genewise_dispersion'])
	return

def makePlotOfPoissonFits(datasets, annotations, outputPath, figSize=(40, 25), nrows=4, ncols=4, title=None, pLimits=(10e-3, 10e4), titleFontSize=16, axesFontSize=14, legendFontSize=10):
	"""
	"""
	
	fig = plt.figure(figsize=figSize)

	outer_grid = fig.add_gridspec(nrows=nrows, ncols=ncols, hspace=0.3, wspace=0.2)

	titleFont={'fontsize' : titleFontSize}
	axesFont={'fontsize' : axesFontSize}

	for i,adata in enumerate(datasets):
	    
		grid_box = outer_grid[i]

		inner_grid = grid_box.subgridspec(2, 1)
		
		# -- Poisson -- 

		ax = fig.add_subplot(inner_grid[0])

		ax.set_title('Poisson', fontdict=titleFont)

		ax.set_xscale('log')
		ax.set_xlim(left=pLimits[0], right=pLimits[1])

		ax.scatter(adata.var['scaled_count_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', label='Observed', rasterized=True);
		ax.scatter(adata.var['scaled_count_mean'],
		       adata.var['poisson_zero_fraction'],
		       ec='w', c='grey', label='Expected', rasterized=True);

		ax.set_ylabel('Fraction zeros', fontdict=axesFont)

		ax.legend(title='Genes', loc='lower left', scatterpoints=3, fontsize=legendFontSize)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Poisson differences

		ax = fig.add_subplot(inner_grid[1])

		ax.set_xscale('log')
		ax.set_xlim(left=pLimits[0], right=pLimits[1])
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['poisson_zero_fraction']
		ax.scatter(adata.var['scaled_count_mean'],
		       difference,
		       c='k', marker='.', label='Genes', rasterized=True)

		ax.set_ylabel('Difference \n(Observed - Expected)', fontdict=axesFont)
		ax.set_xlabel('Mean', fontdict=axesFont)

		ax.legend(loc='lower left', scatterpoints=3, fontsize=legendFontSize)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		# -- Annotation --

		bbox = grid_box.get_position(fig)
		x = (bbox.x0 + bbox.x1) / 2
		y = bbox.y1 + 0.04

		if title:
			figTitle = title
		else: 
			figTitle = adata.uns['name']

		fig.text(x, y, figTitle + '\n' + annotations[i], ha='center', fontdict=titleFont)
	    
	fig.savefig(outputPath, dpi=500, bbox_inches='tight')

def makePlotOfNBFits(datasets, annotations, outputPath, figSize=(40, 25), nrows=4, ncols=4, title=None, limits=(10e-2, 10e3), titleFontSize=16, axesFontSize=14, legendFontSize=10):
	"""
	Makes common and genewise dispersion NB plots in the same manner as the Svensson publication, i.e. no ZI fit.
	"""

	mins = []
	maxs = []

	for adata in datasets:
		difference1 = adata.var['empirical_zero_fraction'] - adata.var['global_zero_fraction']
		difference2 = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
		mins.append(min(difference1.min(), difference2.min()))
		maxs.append(max(difference1.max(), difference2.max()))

	min(mins), max(maxs)
	
	fig = plt.figure(figsize=figSize)

	outer_grid = fig.add_gridspec(nrows=nrows, ncols=ncols, hspace=0.3, wspace=0.2)

	titleFont={'fontsize' : titleFontSize}
	axesFont={'fontsize' : axesFontSize}

	for i,adata in enumerate(datasets):
	    
		grid_box = outer_grid[i]

		inner_grid = grid_box.subgridspec(2, 2)
		# 0 1 2
		# 3 4 5 	
		# -- Global -- 

		ax = fig.add_subplot(inner_grid[0])

		ax.set_title('Common dispersion', fontdict=titleFont)

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])

		ax.scatter(adata.var['empirical_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', label='Observed', rasterized=True);
		ax.scatter(adata.var['empirical_mean'],
		       adata.var['global_zero_fraction'],
		       ec='w', c='grey', label='Expected', rasterized=True);

		ax.set_ylabel('Fraction zeros', fontdict=axesFont)

		ax.legend(title='Genes', loc='lower left', scatterpoints=3, fontsize=legendFontSize)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Global dispersion differences

		ax = fig.add_subplot(inner_grid[2])

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['global_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', label='Genes', rasterized=True)

		ax.set_ylabel('Difference \n(Observed - Expected)', fontdict=axesFont)
		ax.set_xlabel('Mean', fontdict=axesFont)

		ax.legend(loc='lower left', scatterpoints=3, fontsize=legendFontSize)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		# -- Genewise --

		ax = fig.add_subplot(inner_grid[1])

		ax.set_title('Gene-wise dispersion', fontdict=titleFont)

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])

		ax.scatter(adata.var['empirical_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', rasterized=True);
		ax.scatter(adata.var['empirical_mean'],
		       adata.var['genewise_zero_fraction'],
		       ec='w', c='grey', rasterized=True);

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Genewise dispersion differences

		ax = fig.add_subplot(inner_grid[3])

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', rasterized=True)

		ax.set_xlabel('Mean', fontdict=axesFont)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		# -- Annotation --

		bbox = grid_box.get_position(fig)
		x = (bbox.x0 + bbox.x1) / 2
		y = bbox.y1 + 0.05

		if title:
			figTitle = title
		else:
			figTitle = adata.uns['name']

		fig.text(x, y, figTitle + '\n' + annotations[i], ha='center', fontdict=titleFont)
	    
	fig.savefig(outputPath, dpi=500, bbox_inches='tight')

def makePlotOfPoissonAndNB(datasets, annotations, outputPath, figSize=(40, 25), nrows=4, ncols=4, title=None, pLimits=None, nbLimits=None, titleFontSize=16, axesFontSize=14, legendFontSize=10):
	"""
	Combines Poisson, and common and genewise dispersion NB plots.
	"""

	mins = []
	maxs = []

	for adata in datasets:
		difference1 = adata.var['empirical_zero_fraction'] - adata.var['global_zero_fraction']
		difference2 = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
		difference3 = adata.var['empirical_zero_fraction'] - adata.var['poisson_zero_fraction']
		mins.append(min(difference1.min(), difference2.min(), difference3.min()))
		maxs.append(max(difference1.max(), difference2.max(), difference3.max()))

	min(mins), max(maxs)
	
	fig = plt.figure(figsize=figSize)

	outer_grid = fig.add_gridspec(nrows=nrows, ncols=ncols, hspace=0.3, wspace=0.2)

	titleFont={'fontsize' : titleFontSize}
	axesFont={'fontsize' : axesFontSize}

	for i,adata in enumerate(datasets):
	    
		grid_box = outer_grid[i]

		inner_grid = grid_box.subgridspec(2, 3)
		# 0 1 2
		# 3 4 5 	
		
		# -- Poisson -- 

		ax = fig.add_subplot(inner_grid[0])

		ax.set_title('Poisson', fontdict=titleFont)

		ax.set_xscale('log')
		ax.set_xlim(left=pLimits[0], right=pLimits[1])
		# Originally used the scaled_count_mean here, but just complicates presentation
		ax.scatter(adata.var['scaled_count_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', label='Observed', rasterized=True);
		ax.scatter(adata.var['scaled_count_mean'],
		       adata.var['poisson_zero_fraction'],
		       ec='w', c='grey', label='Expected', rasterized=True);

		ax.set_ylabel('Fraction zeros', fontdict=axesFont)

		ax.legend(title='Genes', loc='upper left', scatterpoints=3, fontsize=legendFontSize)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Poisson differences

		ax = fig.add_subplot(inner_grid[3])

		ax.set_xscale('log')
		ax.set_xlim(left=pLimits[0], right=pLimits[1])
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['poisson_zero_fraction']
		ax.scatter(adata.var['scaled_count_mean'],
		       difference,
		       c='k', marker='.', label='Genes', rasterized=True)

		ax.set_ylabel('Difference \n(Observed - Expected)', fontdict=axesFont)
		ax.set_xlabel('Mean', fontdict=axesFont)

		ax.legend(loc='lower left', scatterpoints=3, fontsize=legendFontSize)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		# -- Global -- 

		ax = fig.add_subplot(inner_grid[1])

		ax.set_title('Common dispersion', fontdict=titleFont)

		ax.set_xscale('log')
		ax.set_xlim(left=nbLimits[0], right=nbLimits[1])

		ax.scatter(adata.var['empirical_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', label='Observed', rasterized=True);
		ax.scatter(adata.var['empirical_mean'],
		       adata.var['global_zero_fraction'],
		       ec='w', c='grey', label='Expected', rasterized=True);

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Global dispersion differences

		ax = fig.add_subplot(inner_grid[4])

		ax.set_xscale('log')
		ax.set_xlim(left=nbLimits[0], right=nbLimits[1])
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['global_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', label='Genes', rasterized=True)

		ax.set_xlabel('Mean', fontdict=axesFont)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)


		# -- Genewise --

		ax = fig.add_subplot(inner_grid[2])

		ax.set_title('Gene-wise dispersion', fontdict=titleFont)

		ax.set_xscale('log')
		ax.set_xlim(left=nbLimits[0], right=nbLimits[1])

		ax.scatter(adata.var['empirical_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', rasterized=True);
		ax.scatter(adata.var['empirical_mean'],
		       adata.var['genewise_zero_fraction'],
		       ec='w', c='grey', rasterized=True);

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Genewise dispersion differences

		ax = fig.add_subplot(inner_grid[5])

		ax.set_xscale('log')
		ax.set_xlim(left=nbLimits[0], right=nbLimits[1])
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', rasterized=True)

		ax.set_xlabel('Mean', fontdict=axesFont)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		# -- Annotation --

		bbox = grid_box.get_position(fig)
		x = (bbox.x0 + bbox.x1) / 2
		y = bbox.y1 + 0.04

		if title:
			figTitle = title
		else: 
			figTitle = adata.uns['name']

		fig.text(x, y, figTitle + '\n' + annotations[i], ha='center', fontdict=titleFont)
	    
	fig.savefig(outputPath, dpi=500, bbox_inches='tight')

def makePlotOfGDWithZi(datasets, annotations, outputPath, figSize=(40, 25), nrows=4, ncols=4, title=None, limits=(10e-2,10e4), titleFontSize=16, axesFontSize=14, legendFontSize=10):
	"""
	Makes plots with genewise and ZI fits.
	"""	
	mins = []
	maxs = []

	for adata in datasets:
		difference1 = adata.var['empirical_zero_fraction'] - adata.var['global_zero_fraction']
		difference2 = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
		difference3 = adata.var['empirical_zero_fraction'] - adata.var['genewise_zi_zero_fraction']
		mins.append(min(difference1.min(), difference2.min(), difference3.min()))
		maxs.append(max(difference1.max(), difference2.max(), difference3.max()))

	min(mins), max(maxs)
	
	fig = plt.figure(figsize=figSize)

	outer_grid = fig.add_gridspec(nrows=nrows, ncols=ncols, hspace=0.3, wspace=0.2)

	titleFont={'fontsize' : titleFontSize}
	axesFont={'fontsize' : axesFontSize}

	for i,adata in enumerate(datasets):
	    
		grid_box = outer_grid[i]

		inner_grid = grid_box.subgridspec(2, 2)

		# -- Genewise --

		ax = fig.add_subplot(inner_grid[0])

		ax.set_title('Gene-wise dispersion', fontdict=titleFont)

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])

		ax.scatter(adata.var['empirical_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', rasterized=True);
		ax.scatter(adata.var['empirical_mean'],
		       adata.var['genewise_zero_fraction'],
		       ec='w', c='grey', rasterized=True);

		ax.set_ylabel('Fraction zeros', fontdict=axesFont)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Genewise dispersion differences

		ax = fig.add_subplot(inner_grid[2])

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', rasterized=True)

		ax.set_ylabel('Difference \n(Observed - Expected)', fontdict=axesFont)
		ax.set_xlabel('Mean', fontdict=axesFont)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Zero inflation

		ax = fig.add_subplot(inner_grid[1])

		ax.set_title('ZINB', fontdict=titleFont)

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])

		ax.scatter(adata.var['empirical_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', rasterized=True);
		ax.scatter(adata.var['empirical_mean'],
		       adata.var['genewise_zi_zero_fraction'],
		       ec='w', c='grey', rasterized=True);

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Zero inflation differences

		ax = fig.add_subplot(inner_grid[3])

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['genewise_zi_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', rasterized=True)

		ax.set_xlabel('Mean', fontdict=axesFont)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		# -- Annotation --

		bbox = grid_box.get_position(fig)
		x = (bbox.x0 + bbox.x1) / 2
		y = bbox.y1 + 0.04

		if title:
			figTitle = title
		else: 
			figTitle = adata.uns['name']

		fig.text(x, y, figTitle + '\n' + annotations[i], ha='center', fontdict=titleFont)
	    
	fig.savefig(outputPath, dpi=500, bbox_inches='tight')

def makePlotOfNBWithZi(datasets, annotations, outputPath, figSize=(40, 25), nrows=4, ncols=4, title=None, limits=(10e-2,10e4), titleFontSize=16, axesFontSize=14, legendFontSize=10):

	"""
	Makes plots with genewise ZI fits.
	"""	
	mins = []
	maxs = []

	for adata in datasets:
		difference1 = adata.var['empirical_zero_fraction'] - adata.var['global_zero_fraction']
		difference2 = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
		difference3 = adata.var['empirical_zero_fraction'] - adata.var['genewise_zi_zero_fraction']
		mins.append(min(difference1.min(), difference2.min(), difference3.min()))
		maxs.append(max(difference1.max(), difference2.max(), difference3.max()))

	min(mins), max(maxs)
	
	fig = plt.figure(figsize=figSize)

	outer_grid = fig.add_gridspec(nrows=nrows, ncols=ncols, hspace=0.3, wspace=0.2)

	titleFont={'fontsize' : titleFontSize}
	axesFont={'fontsize' : axesFontSize}

	for i,adata in enumerate(datasets):
	    
		grid_box = outer_grid[i]

		inner_grid = grid_box.subgridspec(2, 3)
		# 0 1 2
		# 3 4 5 	
		# -- Global -- 

		ax = fig.add_subplot(inner_grid[0])

		ax.set_title('Common dispersion', fontdict=titleFont)

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])

		ax.scatter(adata.var['empirical_mean'],
		       adata.var['empirical_zero_fraction'],
		       c='k', label='Observed', rasterized=True);
		ax.scatter(adata.var['empirical_mean'],
		       adata.var['global_zero_fraction'],
		       ec='w', c='grey', label='Expected', rasterized=True);

		ax.set_ylabel('Fraction zeros', fontdict=axesFont)

		ax.legend(title='Genes', loc='lower left', scatterpoints=3, fontsize=8)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Global dispersion differences

		ax = fig.add_subplot(inner_grid[3])

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['global_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', label='Genes', rasterized=True)

		ax.set_ylabel('Difference \n(Observed - Expected)', fontdict=axesFont)
		ax.set_xlabel('Mean', fontdict=axesFont)

		ax.legend(loc='lower left', scatterpoints=3, fontsize=legendFontSize)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)


		# -- Genewise --

		ax = fig.add_subplot(inner_grid[1])

		ax.set_title('Gene-wise dispersion', fontdict=titleFont)

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])

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
		ax.set_xlim(left=limits[0], right=limits[1])
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', rasterized=True)

		ax.set_xlabel('Mean', fontdict=axesFont)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		## Zero inflation

		ax = fig.add_subplot(inner_grid[2])

		ax.set_title('ZINB', fontdict=titleFont)

		ax.set_xscale('log')
		ax.set_xlim(left=limits[0], right=limits[1])

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
		ax.set_xlim(left=limits[0], right=limits[1])
		ax.set_ylim(top=1.0, bottom=-1.0)

		difference = adata.var['empirical_zero_fraction'] - adata.var['genewise_zi_zero_fraction']
		ax.scatter(adata.var['empirical_mean'],
		       difference,
		       c='k', marker='.', rasterized=True)

		ax.set_xlabel('Mean', fontdict=axesFont)

		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)

		# -- Annotation --

		bbox = grid_box.get_position(fig)
		x = (bbox.x0 + bbox.x1) / 2
		y = bbox.y1 + 0.04

		if title:
			figTitle = title
		else: 
			figTitle = adata.uns['name']

		fig.text(x, y, figTitle + '\n' + annotations[i], ha='center', fontdict=titleFont)
	    
	fig.savefig(outputPath, dpi=500, bbox_inches='tight')

def printGOFStatistics(adata):
	print(adata.uns['name'])
	nGenes = adata.var.shape[0]
	zGenes = sum(adata.var.empirical_mean==0)
	dGenes = sum(adata.var.empirical_mean!=0)
	print("Genes: ", nGenes)
	print("Zero genes: ", zGenes)
	print("Non-zero genes: ", dGenes)
	print("Detection rate: %.2f" % (100.0*dGenes/nGenes))
	print(f'{adata.obs.shape[0]} cells, {adata.var.shape[0]} genes')
	difference = adata.var['empirical_zero_fraction'] - adata.var['poisson_zero_fraction']
	absdiff = np.abs(difference)
	print('poisson', (absdiff > 0.01).sum(), (absdiff > 0.05).sum(), (absdiff > 0.1).sum(), (absdiff > 0.2).sum())
	print('poisson %.2f %.2f %.2f %.2f' % (((absdiff > 0.01).sum())*100.0/nGenes, ((absdiff > 0.05).sum())*100.0/nGenes, ((absdiff > 0.1).sum())*100.0/nGenes, ((absdiff > 0.2).sum())*100.0/nGenes))
	difference = adata.var['empirical_zero_fraction'] - adata.var['global_zero_fraction']
	absdiff = np.abs(difference)
	print('common', (absdiff > 0.01).sum(), (absdiff > 0.05).sum(), (absdiff > 0.1).sum(), (absdiff > 0.2).sum())
	print('common %.2f %.2f %.2f %.2f' % (((absdiff > 0.01).sum())*100.0/nGenes, ((absdiff > 0.05).sum())*100.0/nGenes, ((absdiff > 0.1).sum())*100.0/nGenes, ((absdiff > 0.2).sum())*100.0/nGenes))
	difference = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
	absdiff = np.abs(difference)
	print('genewise', (absdiff > 0.01).sum(), (absdiff > 0.05).sum(), (absdiff > 0.1).sum(), (absdiff > 0.2).sum())
	print('genewise %.2f %.2f %.2f %.2f' % (((absdiff > 0.01).sum())*100.0/nGenes, ((absdiff > 0.05).sum())*100.0/nGenes, ((absdiff > 0.1).sum())*100.0/nGenes, ((absdiff > 0.2).sum())*100.0/nGenes))
