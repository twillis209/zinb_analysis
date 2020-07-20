library(zinbwave)

# Writes out (if 
# fig6ad-S13-S14/simAllen_nc100_ratio1_offs-1.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc100_ratio1_offs0.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc100_ratio1_offs2_fitted.rda
# fig6ad-S13-S14/simAllen_nc100_ratio5_offs-1.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc100_ratio5_offs0.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc100_ratio5_offs2_fitted.rda
# fig6ad-S13-S14/simAllen_nc100_ratio10_offs-1.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc100_ratio10_offs0.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc100_ratio10_offs2_fitted.rda
# fig6ad-S13-S14/simAllen_nc1000_ratio1_offs-1.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc1000_ratio1_offs0.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc1000_ratio1_offs2_fitted.rda
# fig6ad-S13-S14/simAllen_nc1000_ratio5_offs-1.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc1000_ratio5_offs0.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc1000_ratio5_offs2_fitted.rda
# fig6ad-S13-S14/simAllen_nc1000_ratio10_offs-1.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc1000_ratio10_offs0.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc1000_ratio10_offs2_fitted.rda
# fig6ad-S13-S14/simAllen_nc10000_ratio1_offs-1.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc10000_ratio1_offs0.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc10000_ratio1_offs2_fitted.rda
# fig6ad-S13-S14/simAllen_nc10000_ratio5_offs-1.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc10000_ratio5_offs0.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc10000_ratio5_offs2_fitted.rda
# fig6ad-S13-S14/simAllen_nc10000_ratio10_offs-1.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc10000_ratio10_offs0.5_fitted.rda
# fig6ad-S13-S14/simAllen_nc10000_ratio10_offs2_fitted.rda
# fig6ad-S13-S14/simZeisel_nc100_ratio1_offs-1.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc100_ratio1_offs0.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc100_ratio1_offs2_fitted.rda
# fig6ad-S13-S14/simZeisel_nc100_ratio5_offs-1.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc100_ratio5_offs0.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc100_ratio5_offs2_fitted.rda
# fig6ad-S13-S14/simZeisel_nc100_ratio10_offs-1.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc100_ratio10_offs0.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc100_ratio10_offs2_fitted.rda
# fig6ad-S13-S14/simZeisel_nc1000_ratio1_offs-1.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc1000_ratio1_offs0.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc1000_ratio1_offs2_fitted.rda
# fig6ad-S13-S14/simZeisel_nc1000_ratio5_offs-1.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc1000_ratio5_offs0.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc1000_ratio5_offs2_fitted.rda
# fig6ad-S13-S14/simZeisel_nc1000_ratio10_offs-1.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc1000_ratio10_offs0.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc1000_ratio10_offs2_fitted.rda
# fig6ad-S13-S14/simZeisel_nc10000_ratio1_offs-1.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc10000_ratio1_offs0.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc10000_ratio1_offs2_fitted.rda
# fig6ad-S13-S14/simZeisel_nc10000_ratio5_offs-1.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc10000_ratio5_offs0.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc10000_ratio5_offs2_fitted.rda
# fig6ad-S13-S14/simZeisel_nc10000_ratio10_offs-1.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc10000_ratio10_offs0.5_fitted.rda
# fig6ad-S13-S14/simZeisel_nc10000_ratio10_offs2_fitted.rda

for (ds in c('Allen', 'Zeisel')) {
	for (nc in c(100, 1000, 10000)) {
		for (b2 in c(1, 5, 10)){
		  for (offs in c(-1.5,0.5,2)){
		    ff = sprintf('fig6ad-S13-S14/sim%s_nc%s_ratio%s_offs%s', ds, nc, b2, offs)
		    load(paste0(ff, '.rda'))
		    fittedSim = lapply(1:4, function(k){
		      lapply(1:length(simData), function(i){
			counts = t(simData[[i]]$counts)
			counts = counts[rowSums(counts) != 0, ]
			ngenes = nrow(counts)
			zinbFit(counts, K = k, commondispersion = FALSE,
				epsilon = ngenes, BPPARAM=SerialParam())
		      })
		    })
		    out = paste0(ff, '_fitted.rda')
		    save(fittedSim, file = out)
		  }
		}
	}
}
