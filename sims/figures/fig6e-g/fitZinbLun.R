library(zinbwave)

# Writes out 
# fig6e-g/simLun_100_ziadd0_fitted.rda
# fig6e-g/simLun_100_ziadd0.33_fitted.rda
# fig6e-g/simLun_100_ziadd0.67_fitted.rda
# fig6e-g/simLun_1000_ziadd0_fitted.rda
# fig6e-g/simLun_1000_ziadd0.33_fitted.rda
# fig6e-g/simLun_1000_ziadd0.67_fitted.rda
# fig6e-g/simLun_10000_ziadd0_fitted.rda
# fig6e-g/simLun_10000_ziadd0.33_fitted.rda
# fig6e-g/simLun_10000_ziadd0.67_fitted.rda

fitAllSimLun<-function() {
	for (ncells in c(100, 1000, 10000)){
	  for (add in c('_ziadd0', '_ziadd0.33', '_ziadd0.67')){
		fitSimData(ncells, add)
	  }
	}
}

fitSimData<-function(ncells, add) {
    fileName = sprintf('simLun_%s%s', ncells, add)
    load(paste0(fileName, '.rda'))
    fittedSim = lapply(1:length(simData), function(i){
      counts = simData[[i]]$counts
      counts = counts[rowSums(counts) != 0, ]
      zinbFit(counts, K = 2, commondispersion = FALSE)
    })
    out = paste0(fileName, '_fitted.rda')
    save(fittedSim, file = out)
}

fitSimData(10000, '_ziadd0.33')
