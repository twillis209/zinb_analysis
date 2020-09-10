library(digest)
library(zinbwave)

# Defaulting to common dispersion here as that is true of zinbFit
fitZinbAndCache<-function(countMatrix, cacheDirPath="zinbCache", K=2, epsilon=1000, commonDispersion=T) {
  # fit zinb (if you already fitted zinb, it is cached)
  d = digest::digest(countMatrix, "md5")
 
  cachePath = paste0(cacheDirPath, '/', d)
 
  fileZinb = sprintf("%s_%d_%d_%d_zinb.rda", cachePath, K, epsilon, commonDispersion)
 
  if (!file.exists(fileZinb)){
    print('run ZINB')
 
    zinb <- zinbFit(countMatrix, K = K, commondispersion = commonDispersion, epsilon = epsilon)
    save(zinb, file = fileZinb)
  }else{
    load(fileZinb)
  }
  return(zinb)
}

# Different signature for design matrix
fitZinbAndCache<-function(countMatrix, X, cacheDirPath="zinbCache", K=2, epsilon=1000, commonDispersion=T) {
  # fit zinb (if you already fitted zinb, it is cached)
  d = digest(countMatrix, "md5")
 
  cachePath = paste0(cacheDirPath, '/', d)
 
  fileZinb = sprintf("%s_%d_%d_%d_zinb.rda", cachePath, K, epsilon, commonDispersion)
 
  if (!file.exists(fileZinb)){
    print('run ZINB')
 
    zinb <- zinbFit(countMatrix, X=X, K = K, commondispersion = commonDispersion, epsilon = epsilon)
    save(zinb, file = fileZinb)
  }else{
    load(fileZinb)
  }
  return(zinb)
}

# Different signature for design matrices
fitZinbAndCache<-function(countMatrix, X, V, cacheDirPath="zinbCache", K=2, epsilon=1000, commonDispersion=T) {
  # fit zinb (if you already fitted zinb, it is cached)
  d = digest(countMatrix, "md5")
 
  cachePath = paste0(cacheDirPath, '/', d)
 
  fileZinb = sprintf("%s_%d_%d_%d_zinb.rda", cachePath, K, epsilon, commonDispersion)
 
  if (!file.exists(fileZinb)){
    print('run ZINB')
 
    zinb <- zinbFit(countMatrix, X=X, V=V, K = K, commondispersion = commonDispersion, epsilon = epsilon)
    save(zinb, file = fileZinb)
  }else{
    load(fileZinb)
  }
  return(zinb)
}

wrapRzifa <- function(Y, block = TRUE, k=2, cacheDirPath="zifaCache"){
  # wrapper R function for ZIFA.
  # md5 hashing and temporary files are used not to re-run zifa 
  # if it has already be run on this computer.
  d = digest(Y, "md5")

  cachePath = paste0(cacheDirPath, '/', d)

  if (!file.exists(paste0(cachePath, "_", k, '_zifa.csv'))){
    write.csv(Y, paste0(cachePath, ".csv"))
    print('run ZIFA')
    bb = ifelse(block, '-b ', '')
    cmd = sprintf('python run_zifa.py -d %d %s %s.csv %s_%d_zifa.csv', k, bb, cachePath, cachePath, k)
    system(cmd)
  }
  read.csv(sprintf("%s_%d_zifa.csv", cachePath, k), header=FALSE)
}
