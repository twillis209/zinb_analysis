library(zinbwave)

# Writes out cpuTime.rda

simZeisel_nc%s_ratio1_offs2
cpuTime = lapply(c(50, 100, 500, 1000, 5000, 10000), function(nc){
  fileName = sprintf('simZeisel_nc%s_ratio1_offs2', nc)
  load(paste0(fileName, ".rda"))
  tt = lapply(1:10, function(j){
    counts = t(simData[[j]]$counts)
    counts = counts[rowSums(counts) > 5, colSums(counts) > 5 ] 
    system.time(zinbFit(counts, K = 2, commondispersion = FALSE,
                        epsilon = 1000))
  })
  tt
})
save(cpuTime, file = 'cpuTime.rda')
