library(BiocParallel)
library(zinbwave)

# Writes out (depending on CL args)
# fig6ad-S13-S14/sim(ds)_nc(nc)_ratio(b2)_offs(offs)_k1.rda
# fig6ad-S13-S14/sim(ds)_nc(nc)_ratio(b2)_offs(offs)_k2.rda
# fig6ad-S13-S14/sim(ds)_nc(nc)_ratio(b2)_offs(offs)_k3.rda
# fig6ad-S13-S14/sim(ds)_nc(nc)_ratio(b2)_offs(offs)_k4.rda
# fig6ad-S13-S14/sim(ds)_nc(nc)_ratio(b2)_offs(offs)_fitted.rda

fitSim<-function(ds, nc, b2, offs, BPPARAM=BiocParallel::bpparam()) {
	ff = sprintf('sim%s_nc%s_ratio%s_offs%s', ds, nc, b2, offs)
	load(paste0(ff, '.rda'))
	fittedSim = lapply(1:4, function(k){
	  print(k)
	  aa = lapply(1:2, function(i){
	    print(i)
	    set.seed(73839)
	    counts = t(simData[[i]]$counts)
	    counts = counts[rowSums(counts) != 0, colSums(counts) != 0]
	    ngenes = nrow(counts)
	    zinbFit(counts, K = k, commondispersion = FALSE, verbose = T, epsilon = ngenes, BPPARAM=SerialParam())
	  })
	  save(aa, file = sprintf('%s_k%s.rda', ff, k))
	  aa
	})
	out = paste0(ff, '_fitted.rda')
	save(fittedSim, file = out)
}
