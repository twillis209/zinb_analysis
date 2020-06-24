library(scRNAseq)
library(zinbwave)
library(BiocParallel)
library(SingleCellExperiment)
library(magrittr)
library(tictoc)
# Just a quick script to see if I am using the BiocParallel interface correctly.

# Following steps in zinbwave vignette

#register(MulticoreParam(2))

allenSCE<-ReprocessedAllenData()

hvgFilter<-rowSums(assay(allenSCE)>5)>5

allenSCE_hvgs<-allenSCE[hvgFilter,]

# Perhaps not available for this version of SingleCellExperiment?
#rowSubset(allenSCE, "hvgs")<-rowSums(assay(allenSCE)>5)>5

assay(allenSCE_hvgs) %>% log1p %>% rowVars-> vars
names(vars)<-rownames(allenSCE_hvgs)
vars<-sort(vars, decreasing=T)

allenSCE_hvgs<-allenSCE_hvgs[names(vars)[1:100],]

timeZinbwave<-function(){
	tic()
	allenSCE_zinb<-zinbwave(allenSCE_hvgs, K=2, epsilon=100, which_assay="tophat_counts")
	toc(log=T, quiet=T)
}
suppressWarnings({
replicate(5, timeZinbwave())
tic.log()
print(unlist(tic.log()))
tic.clearlog()
})
