library(BiocParallel)
library(zinbwave)

# Writes out cpuTime.rda

cpuTime = bplapply(c(50, 100, 500, 1000, 5000, 10000), function(nc){
  fileName = sprintf('simZeisel_nc%s_ratio1_offs2', nc)
  load(paste0(fileName, ".rda"))
  tt = lapply(1:10, function(j){
    counts = t(simData[[j]]$counts)
    counts = counts[rowSums(counts) > 5, colSums(counts) > 5 ] 
    system.time(zinbFit(counts, K = 2, commondispersion = FALSE,
                        epsilon = 1000, BPPARAM=SerialParam()))
  })
  tt
}, BPPARAM=SerialParam())
save(cpuTime, file = 'cpuTime.rda')
