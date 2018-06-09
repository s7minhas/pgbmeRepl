# workspace ###############################
setwd('~/Research/pgbmeRepl/main')
rm(list=ls())
library(tidyr)
library(reshape2)
library(sbgcop)
library(doParallel)
library(foreach)

# helpers
char = function(x){as.character(x)}
num = function(x){as.numeric(char(x))}
cntr <- function(x) (x - mean(c(x), na.rm = TRUE))/sd(c(x), na.rm = TRUE)
trim = function (x) { gsub("^\\s+|\\s+$", "", x) }
################################

# load data files and assemble for imp ###############################
load('rawData.rda') # bit.acc.t, covData, dist.norm

# only covData needs to be imputed
covDataL = do.call('rbind', covData)
forSbg = covDataL[,-c(1:3)]
sbgMod = sbgcop.mcmc(forSbg, nsamp=1000, odens=1, seed=6886)

# pull out a few imp datasets
set.seed(6886)
impIters = sample(501:1000, 11)
sbgImps = lapply(impIters, function(i){
	x = sbgMod$Y.impute[,,i] ; colnames(x) = colnames(sbgMod$Y.pmean)
	return( data.frame(covDataL[,1:3], x, row.names=NULL) ) })
sbgImps[[1]] = data.frame(covDataL[,1:3], sbgMod$Y.pmean,row.names=NULL)
################################

# design structs ###############################
cl=makeCluster(2) ; registerDoParallel(cl)
shhh <- foreach(t = c('1995','2010')) %dopar% {

	# iterate over every imp
	xData = lapply(sbgImps, function(x){
		xN = x[x$year==num(t),] ; xN$gdpCapLog = log(xN$gdpCap + 1)
		rownames(xN) = xN$Country.Code ; xN = xN[rownames(bit.acc.t[[t]]),]

		# dyad covars
		xDyad = array(NA, dim=c(dim(bit.acc.t[[t]]), 5), 
			dimnames=list(rownames(bit.acc.t[[t]]), rownames(bit.acc.t[[t]]), 
				c('udsDiff', 'lawOrderDiff', 'gdpCapLogDiff', 'oecdJoint', 'dist')))
		for(i in rownames(xDyad)){
			xDyad[i,,'udsDiff'] = abs( xN[i,'uds_median'] - xN[,'uds_median'])
			xDyad[i,,'lawOrderDiff'] = abs( xN[i,'lawOrder'] - xN[,'lawOrder'])
			xDyad[i,,'gdpCapLogDiff'] = abs( xN[i,'gdpCapLog'] - xN[,'gdpCapLog'])
			xDyad[i,,'oecdJoint'] = xN[i,'oecd'] * xN[,'oecd'] }
		xDyad[,,'dist'] = dist.norm[[t]]
		for(p in 1:dim(xDyad)[3]){ xDyad[,,p] = cntr(xDyad[,,p]) }

		# nodal covars
		xN = data.matrix(xN[,c('fdiGdp','dispC','gdpCapGr','ptaCumul')])
		xN = apply(xN, 2, cntr)

		# red
		return(list(xDyad=xDyad,xNode=xN)) })

	# save
	save(bit.acc.t, xData, 
		file=paste0('modelData',t,'.rda'))
}
stopCluster(cl)
################################