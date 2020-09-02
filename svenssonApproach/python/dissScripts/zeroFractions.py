import anndata


datasets = [anndata.read(x) for x in ['../../Data/output/sheng/noFilterSetOne293T_sAndp_exonU.h5ad', '../../Data/output/sheng/noFilter293T_singleCell_exonU.h5ad']]

for x in datasets:
        print(x.uns['name'])
        print("Genes: ", x.var.shape[0])
        print(sum(x.var.empirical_mean==0))
