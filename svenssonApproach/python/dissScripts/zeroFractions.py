import anndata


datasets = [anndata.read(x) for x in ['../../Data/output/sheng/ercc.h5ad', '../../Data/output/sheng/noFilterSetOne293T_sAndp_exonU.h5ad', '../../Data/output/sheng/noFilter293T_singleCell_exonU.h5ad']]
types = ["ercc", "sAndp", "singleCell"]
for i,x in enumerate(datasets):
	print(types[i])
	print(x.uns['name'])
	nGenes = x.var.shape[0]
	zGenes = sum(x.var.empirical_mean==0)
	dGenes = sum(x.var.empirical_mean!=0)
	print("Genes: ", x.var.shape[0])
	print("Zero genes: ", sum(x.var.empirical_mean==0))
	print("Non-zero genes: ", sum(x.var.empirical_mean != 0))
	print("Detection rate: %.2f" % (100.0*dGenes/nGenes))
