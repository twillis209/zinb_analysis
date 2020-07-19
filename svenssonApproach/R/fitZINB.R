library(countreg)

kolodDAF<-read.table('../Data/kolodziejczyk_counttable_es.csv', sep=' ', header=T)

erccDAF<-kolodDAF[grepl('^ERCC', rownames(kolodDAF)),]
# ercc2 is ERCC-00002
ercc2<-t(erccDAF[1,])
#init.theta=ifelse(var(ercc2) > mean(ercc2), (mean(ercc2)^2)/(var(ercc2)-mean(ercc2)), 10)[1]
nbModel<-glm.nb(ercc2 ~ 1)
mu=exp(nbModel$coefficients[1])

# 'size' in Svensson's code is the r parameter, which is the reciprocal of the dispersion parameter
# theta appears to be the reciprocal of alpha in the 'direct' parameterisation
# The initial 'size' estimate employed by Svensson is the reciprocal of the initial theta estimate used by glm.nb, but theta ultimately converges on a different value here

zinbModel<-zeroinfl(ercc2 ~ 1|1, dist="negbin")
# For ERCC-00002, zeroinfl produced the same theta estimate for the NB component as glm.nb
erccZinbModels<-list()
for(i in seq(1:nrow(erccDAF))) {
	print(rownames(erccDAF)[i])
	erccZinbModels[[rownames(erccDAF)[i]]]<-tryCatch({zeroinfl(t(erccDAF[i,])~1|1, dist="negbin")}, warning=function(w) { NULL}, error=function(e) { NULL})
}

subErccDAF<-erccDAF[names(erccZinbModels),]

erccStatistics=data.frame(empirical_mean=apply(subErccDAF, 1, mean),
			empirical_variance=apply(subErccDAF, 1, var),
			empirical_zero_fraction=apply(subErccDAF, 1, function(x) sum(x==0)/length(x)),
			ml_mean=sapply(erccZinbModels, function(x) exp(x$coefficients$count)),
			genewise_dispersion=sapply(erccZinbModels, function(x) 1/x$theta),
			genewise_zi_zero_fraction=sapply(erccZinbModels, function(x) dzinbinom(0, mu=exp(x$coefficients$count), 
											theta=x$theta, pi=(exp(x$coef$zero)/(1+exp(x$coef$zero))))))

write.csv(erccStatistics, file="../Data/output/kolodERCCzi.csv")
