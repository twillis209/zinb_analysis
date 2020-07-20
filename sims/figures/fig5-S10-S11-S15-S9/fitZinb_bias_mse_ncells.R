library(zinbwave)
library(BiocParallel)

# Writes out: 
# 
# simZeisel_nc10000_ratio1_offs2_fitted.rda
# simZeisel_nc1000_ratio1_offs2_fitted.rda
# simZeisel_nc100_ratio1_offs2_fitted.rda
# simZeisel_nc5000_ratio1_offs2_fitted.rda
# simZeisel_nc500_ratio1_offs2_fitted.rda
# simZeisel_nc50_ratio1_offs2_fitted.rda

fitSim<-function(nc, ds, b2, offs, eps, BPPARAM=BiocParallel::bpparam()) {
  pp = sprintf('simZeisel_nc%s_ratio1_offs2', nc)
  load(paste0(pp, ".rda"))
  fittedSim = lapply(1:length(simData), function(j){
    counts = t(simData[[j]]$counts)
    counts = counts[rowSums(counts) != 0, colSums(counts) != 0]
    zinbFit(counts, K = 2, commondispersion = FALSE, epsilon = eps, BPPARAM=BPPARAM)
  })
  out = paste0(pp, '_fitted.rda')
  save(fittedSim, file = out)
}

nc = c(50, 100, 500, 1000, 5000, 10000)
ds = 'Zeisel'
b2 = 1
offs = 2
eps = 1000
lapply(nc, function(x) fitSim(nc=x, ds=ds, b2=b2, offs=offs, eps=eps, BPPARAM=SerialParam()))
