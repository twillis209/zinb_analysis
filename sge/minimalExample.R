library(BiocParallel)

## Pi approximation
piApprox <- function(n) {
nums <- matrix(runif(2 * n), ncol = 2)
d <- sqrt(nums[, 1]^2 + nums[, 2]^2)
4 * mean(d <= 1)
}

#param <- BatchtoolsParam(workers=1, cluster="sge", template=paste(Sys.getenv("HOME"), "sge-simple.tmpl", sep="/"), progressbar=T,
#		resources=list(job.name="minimalExample_child",
#				h_rt="0:05:00",
#				mem_free="8G",
#				h_vmem="8G",
#				smp="8"),
#		registryargs=batchtoolsRegistryargs(
#				file.dir=paste("/nobackup", Sys.getenv("USER"), "minimalExample", sep="/"),
#				work.dir=paste(Sys.getenv("HOME"), "zinb_analysis/sims/figures/fig5-S10-S11-S15-S9", sep="/")	
#				))
#logdir=paste(Sys.getenv("HOME"), "logs", sep="/")

param<-SnowParam(workers=16, progressbar=T, log=T, jobname="minimalExample")
bpstart(param)
## Run parallel job
result <- bplapply(rep(10e5, 1e4), piApprox, BPPARAM=param)
bpstop(param)
print("Done")
