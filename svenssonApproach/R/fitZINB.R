library(pscl)

kolodDAF<-read.table('../Data/kolodziejczyk_counttable_es.csv', sep=' ', header=T)

erccDAF<-kolodDAF[grepl('^ERCC', rownames(kolodDAF)),]

counts<-data.frame(counts=t(erccDAF)[, "ERCC-00002"])

nbMod<-zeroinfl(counts ~ 1, data=counts, link="logit", dist="negbin")
