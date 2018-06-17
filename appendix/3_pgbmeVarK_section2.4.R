# workspace ###############################
rm(list=ls())
path <- '/home/minhas/' # ubuntu path format for ec2
# path <- '~/Research/pgbmeRepl/' # example path format for mac
aPath <- paste0(path, 'appendix/')
mPath <- paste0(path, 'main/')
setwd(aPath)

# install packages
toInstall <- c(
	'magic', 'msm', 'lme4', 'mnormt', 'abind', 
	'foreach', 'doParallel', 'PRROC'
	)
for(pkg in toInstall){
  if(!pkg %in% installed.packages()[,1]){
    install.packages(pkg) } }

# load libraries
library(magic)
library(msm)
library(lme4)
library(mnormt)
library(abind)
library(foreach)
library(doParallel)
library(PRROC)

# helpers
library(pgbme)
char = function(x){as.character(x)}
num = function(x){as.numeric(char(x))}
cntr <- function(x) (x - mean(c(x), na.rm = TRUE))/sd(c(x), na.rm = TRUE)
trim = function (x) { gsub("^\\s+|\\s+$", "", x) }
################################

# get data ###############################
load( paste0(mPath, 'modelData2012.rda') )

# sample from posterior of imputed datasets
y <- data.matrix(bit.acc.t[[ '2012' ]])
cntryKey = data.frame(cbind(name=rownames(y), id=1:nrow(y)))
rownames(y) = colnames(y) = cntryKey$id
y <- apply(mat.vect(y), 1, prod)

#
xDyad = xData[[1]]$xDyad
for(p in 1:dim(xDyad)[3]){ diag(xDyad[,,p])=NA }
dimnames(xDyad)[[1]] = dimnames(xDyad)[[2]] = cntryKey$id

#
xNode = xData[[1]]$xNode
rownames(xNode) = cntryKey$id
rm(xData)
################################

# run pgbme varying k ###############################
if(!all(paste0('results2012_k', 1:3, '.rda') %in% list.files())){
	cores = 3
	cl=makeCluster(cores) ; registerDoParallel(cl)
	shhh <- foreach(rank = 1:3, 
		.packages=c('lme4','magic','msm','mnormt','pgbme')
		) %dopar% {
		set.seed(6886)
		est <- pgbme(
			y = y, Xd = xDyad, Xs = xNode, Xr = xNode, 
			k = rank, rho.calc = FALSE,
			NS = 2e+4, burn = 1e+4, odens = 10, 
			xInclImpList=FALSE
			)
		save(est, file=paste0('results2012_k',rank,'.rda'))
	}
	stopCluster(cl)
}
################################

# pull out yhats ###############################
predDF <- lapply(1:3, function(rank){
	load(paste0('results2012_k',rank,'.rda'))
	yhat <- calc_yhat(est, FALSE); rm(est)
	yhat <- apply(pnorm(yhat), 1, mean)
	yhat <- matrix(yhat, sqrt(length(yhat)), sqrt(length(yhat)))
	diag(yhat) = NA ; pgbmeProb = c(yhat) ; rm(yhat)
	return(pgbmeProb) })
predDF <- data.frame(do.call('cbind', predDF))
names(predDF) = paste0('k_',1:3)
################################

# calc auc stats ###############################
predDF$y = c(data.matrix(bit.acc.t[[ '2012' ]]))
predDF = na.omit(predDF)

# run perf checks
perfStats = lapply(names(predDF)[1:3], function(k){
	aucROC=roc.curve(predDF[predDF$y==1,k], predDF[predDF$y==0,k])$auc
	aucPR=pr.curve(predDF[predDF$y==1,k], predDF[predDF$y==0,k])$auc.integral
	return( c('AUC (ROC)'=aucROC,'AUC (PR)'=aucPR) )
})
perfStats = do.call('rbind', perfStats)
rownames(perfStats) = names(predDF)[1:3]
sink(file='tableA3.txt')
print(round(perfStats, 2))
sink()
################################