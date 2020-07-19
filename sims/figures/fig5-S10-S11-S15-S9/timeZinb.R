library(BiocParallel)
library(zinbwave)

#param <- BatchtoolsParam(workers=1, cluster="sge", template=paste(Sys.getenv("HOME"), "sge-simple.tmpl", sep="/"),
#		resources=list(job.name="timeZinb_child",
#				h_rt="1:00:00",
#				mem_free="12G",
#				h_vmem="12G",
#				smp="8"),
#		registryargs=batchtoolsRegistryargs(
#				file.dir=paste("/nobackup", Sys.getenv("USER"), "timeZinb", sep="/"),
#				work.dir=paste(Sys.getenv("HOME"), "zinb_analysis/sims/figures/fig5-S10-S11-S15-S9", sep="/")	
#				))

# Writes out cpuTime.rda

param<-MulticoreParam(workers=6, progressbar=T, log=T, jobname="timeZinb")

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
}, BPPARAM=param)
save(cpuTime, file = 'cpuTime.rda')
