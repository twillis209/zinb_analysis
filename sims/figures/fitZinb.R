library(scRNAseq)
library(zinbwave)
library(mclust)
library(digest)
library(RColorBrewer)
library(MASS)
library(magrittr)
library(BiocParallel)

load("../../real_data/allen/allen.rda")
cols = brewer.pal(8, "Set1")
prefilter = allen[grep("^ERCC-", rownames(allen), invert = TRUE),
		   which(colData(allen)$Core.Type=="Core")]
filterGenes = apply(assay(prefilter) > 5, 1, sum) >= 5
postfilter = prefilter[filterGenes, ]
bioIni =  as.factor(colData(postfilter)$driver_1_s)
core = assay(postfilter)
colIni = cols[bioIni]

param<-MulticoreParam(workers=6, log=T)
zinb <- zinbFit(core, K = 2, commondispersion = F, epsilon = 1000, BPPARAM=param)
