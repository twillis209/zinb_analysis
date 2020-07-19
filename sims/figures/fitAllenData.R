source('simFunction.R')

#param <- BatchtoolsParam(workers=8, 
#			cluster="sge", 
#			template=paste(Sys.getenv("HOME"), "sge-simple.tmpl", sep="/"),
#			progressbar=T,
#			log=T,
#			resources=list(job.name="fitAllenData_child",
#				h_rt="0:05:00",
#				mem_free="4G",
#				h_vmem="4G",
#				smp="1"),
#			registryargs=batchtoolsRegistryargs(
#				file.dir=paste("/nobackup", Sys.getenv("USER"), "fitAllenData", sep="/"),
#				work.dir=paste(Sys.getenv("HOME"), "zinb_analysis/sims/figures", sep="/")	
#				))

#param <- BatchtoolsParam(workers=4, 
#			cluster="multicore", 
#			progressbar=T,
#			log=T,
#			registryargs=batchtoolsRegistryargs(
#				file.dir=paste("/nobackup", Sys.getenv("USER"), "simulateAllenData", sep="/"),
#				work.dir=paste(Sys.getenv("HOME"), "zinb_analysis/sims/figures", sep="/")	
#				))

param<-MulticoreParam(workers=4, tasks=4, log=T, progressbar=T)

loadAndFilterAllenData()

fitZinbAndCache(core, cacheDirPath="zinbCache4")
