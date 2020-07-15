from glob import glob
import anndata
import pandas as pd
import matplotlib.pyplot as plt

datasets = [
    anndata.read('../Data/output/allen.h5ad'),
    anndata.read('../Data/output/fletcher.h5ad'),
    anndata.read('../Data/output/kolod.h5ad'),
    anndata.read('../Data/output/zheng.h5ad'),
]

annotation = {
    'Tasic et al. 2016': 'ERCC spike-ins',
    'Fletcher et al. 2017': 'ERCC spike-ins',
    'Kolodziejczyk et al. 2015': 'ERCC spike-ins',
    'Zheng et al. 2017': 'ERCC spike-ins'
}

location = {
    'Tasic et al. 2016': 0,
    'Fletcher et al. 2017': 1,
    'Kolodziejczyk et al. 2015': 2 ,
    'Zheng et al. 2017': 3 
}

mins = []
maxs = []
for adata in datasets:
    difference1 = adata.var['empirical_zero_fraction'] - adata.var['global_zero_fraction']
    difference2 = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
    mins.append(min(difference1.min(), difference2.min()))
    maxs.append(max(difference1.max(), difference2.max()))

min(mins), max(maxs)

fig = plt.figure(figsize=(15, 25))

outer_grid = fig.add_gridspec(4, 2, hspace=0.4, wspace=0.3)

for adata in datasets:
    
    i = location[adata.uns['name']]
    grid_box = outer_grid[i]
    
    inner_grid = grid_box.subgridspec(2, 2)
    
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
    
    ## _
    
    ax = fig.add_subplot(inner_grid[2])
    
    ax.set_xscale('log')
    ax.set_xlim(left=2e-4, right=1e4)
    ax.set_ylim(top=0.9, bottom=-0.5)
    
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
    
    ## _
    
    ax = fig.add_subplot(inner_grid[3])
    
    ax.set_xscale('log')
    ax.set_xlim(left=2e-4, right=1e4)
    ax.set_ylim(top=0.9, bottom=-0.5)
    
    difference = adata.var['empirical_zero_fraction'] - adata.var['genewise_zero_fraction']
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
    
    fig.text(x, y, adata.uns['name'] + '\n' + annotation[adata.uns['name']], ha='center', fontsize=12)
    
fig.savefig('../Figures/nb_plot.pdf', dpi=500, bbox_inches='tight')

