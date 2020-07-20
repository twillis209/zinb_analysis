library(countreg)

kolodDAF<-read.table('../Data/kolodziejczyk_counttable_es.csv', sep=' ', header=T)

prefixes=c('ola_mES_lif', 'ola_mES_2i', 'ola_mES_a2i', 'ola_mES_lif_1', 'ola_mES_lif_2', 'ola_mES_lif_3', 'ola_mES_2i_2', 'ola_mES_2i_3', 'ola_mES_2i_4','ola_mES_2i_5', 'ola_mES_a2i_2', 'ola_mES_a2i_3')

kolodDFs<-lapply(prefixes, function(x) kolodDAF[!grepl('^ERCC|^_', rownames(kolodDAF)),grepl(x, names(kolodDAF))])
kolodDFs<-lapply(kolodDFs, function(x) x[apply(x, 1, function(y) any(y==0)),])
names(kolodDFs)<-prefixes

# Can't fit ZINB models to features with no zero counts; how many are like this?

# 'size' in Svensson's code is the r parameter, which is the reciprocal of the dispersion parameter
# theta appears to be the reciprocal of alpha in the 'direct' parameterisation
# The initial 'size' estimate employed by Svensson is the reciprocal of the initial theta estimate used by glm.nb, but theta ultimately converges on a different value here

zinbModels<-list()
for(j in 1:length(kolodDFs)) {
	zinbModels[[names(kolodDFs)[[j]]]]<-list()

	for(i in 1:nrow(kolodDFs[[j]])) {
#	for(i in 1:10) {
		print(rownames(kolodDFs[[j]])[i])
	zinbModels[[names(kolodDFs)[[j]]]][[rownames(kolodDFs[[j]])[i]]]<-tryCatch({zeroinfl(t(kolodDFs[[j]][i,])~1|1, dist="negbin")}, warning=function(w) { message(w); NULL}, error=function(e) { message(e); NULL})
	}
}

# Only keep those genes for which we could fit models
modelledKolodDFs<-list()
for(i in 1:length(kolodDFs)) {
	modelledKolodDFs[[names(kolodDFs)[[i]]]]<-kolodDFs[[names(kolodDFs)[[i]]]][names(zinbModels[[names(kolodDFs)[[i]]]]),]
}

statsDFs<-list()

for(i in 1:length(modelledKolodDFs)) {
	mkdf<-modelledKolodDFs[[i]]
	name<-names(modelledKolodDFs)[[i]]
	statsDFs[[name]]<-data.frame(empirical_mean=apply(mkdf, 1, mean),
			empirical_variance=apply(mkdf, 1, var),
			empirical_zero_fraction=apply(mkdf, 1, function(x) sum(x==0)/length(x)),
			ml_mean=sapply(zinbModels[[name]], function(x) exp(x$coefficients$count)),
			genewise_dispersion=sapply(zinbModels[[name]], function(x) 1/x$theta),
			genewise_zi_zero_fraction=sapply(zinbModels[[name]], function(x) dzinbinom(0, mu=exp(x$coefficients$count), 
											theta=x$theta, pi=(exp(x$coef$zero)/(1+exp(x$coef$zero))))))
	write.csv(statsDFs[[name]], file=paste("../Data/output/kolodZi", paste0(name, ".csv"), sep="_"))
}
