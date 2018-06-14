# workspace ###############################
rm(list=ls())
path <- '~/Research/pgbmeRepl/main/'
setwd(path)

# load libraries
library(magic)
library(msm)
library(lme4)
library(mnormt)
library(abind)
library(reshape2)
library(dplyr)
library(tidyr)
library(PRROC)
library(foreach)
library(doParallel)

# helpers
source("pgbme.R")
char = function(x){as.character(x)}
num = function(x){as.numeric(char(x))}
cntr <- function(x) (x - mean(c(x), na.rm = TRUE))/sd(c(x), na.rm = TRUE)
trim = function (x) { gsub("^\\s+|\\s+$", "", x) }
################################

# get data ###############################
load('modelData2012.rda')
y <- data.matrix(bit.acc.t[[ '2012' ]])
cntryKey = data.frame(cbind(name=rownames(y), id=1:nrow(y)))
rownames(y) = colnames(y) = cntryKey$id

#
xDyad = xData[[1]]$xDyad
for(p in 1:dim(xDyad)[3]){ diag(xDyad[,,p])=NA }
dimnames(xDyad)[[1]] = dimnames(xDyad)[[2]] = cntryKey$id

#
xNode = xData[[1]]$xNode
rownames(xNode) = cntryKey$id
rm(xData)

# create na mat
set.seed(6886)
nFolds = 5
foldMat = matrix(
	sample(1:nFolds, 
		nrow(bit.acc.t[[ '2012' ]])^2, replace=TRUE),
	nrow=nrow(y))
################################

# run pgbme ###############################
fName <- paste0('resultsPGBME_outPerf_folds.rda')
if(
	!fName %in% list.files()
	){
	cores = 5
	cl=makeCluster(cores) ; registerDoParallel(cl)
	pgbmeResults <- foreach(
		fold = 1:nFolds, 
		.packages=c('lme4','magic','msm','mnormt')) %dopar% {

		#
		source('pgbme.R')

		# add NAs to dv
		naMat = foldMat ; naMat[naMat==fold] <- NA ;
		naMat[naMat!=fold] <- 1 ; diag(naMat) = NA ;
		yMiss = y * naMat
		yMiss <- apply(mat.vect(yMiss), 1, prod)

		set.seed(6886)
		pgbmeEst <- pgbme(
			y=yMiss, Xd=xDyad, Xs=xNode, Xr=xNode, 
			k = 2, rho.calc = FALSE,
			NS = 2e+4, burn = 1e+4, odens = 10,
			xInclImpList=FALSE, seed=6886
			)
		return(pgbmeEst)
	}
	stopCluster(cl)
	save(pgbmeResults, file=fName)
}
################################

# run gbme ###############################
if(
	!all(file.exists(
		paste0(
			'gbme_outPerf/fold',1:nFolds,'/OUT'
			)))
	){
	cores = 5
	cl=makeCluster(cores) ; registerDoParallel(cl)
	shh<-foreach(
		fold = 1:nFolds, 
		.packages=c('lme4')) %dopar% {

		#
		source('gbme.R')

		#
		naMat = foldMat ; naMat[naMat==fold] <- NA ;
		naMat[naMat!=fold] <- 1 ; diag(naMat) = NA ;
		yMiss = y * naMat

		#
		setwd(paste0(path, 'gbme_outPerf/fold',fold,'/'))
		set.seed(6886)
		est <- gbme(
			Y = yMiss, Xd = xDyad, Xs = xNode, 
			fam='binomial', directed=FALSE, k = 2, 
			N = matrix(1,nrow(y),nrow(y)),
			NS = 2e+4, odens = 10, seed=6886,	
			owrite=TRUE, awrite=TRUE
			)
	}
	stopCluster(cl)
}
################################

# org results, run glm, calc auc stats   ###############################
load(paste0(path, 'resultsPGBME_outPerf_folds.rda'))
predList = lapply(1:nFolds, function(fold){

	# read in gbme model results
	source('gbme.R')
	ddesign <- calc_yhat_gbme(y, xDyad, xNode, 
		paste0(path, 'gbme_outPerf/fold',fold,'/'))

	# pull in pgbme results
	source("pgbme.R")
	pgbmeEst <- pgbmeResults[[fold]]
	yhat <- calc_yhat(pgbmeEst, FALSE); rm(pgbmeEst)
	yhat <- apply(pnorm(yhat), 1, mean)
	yhat <- matrix(yhat, sqrt(length(yhat)), sqrt(length(yhat)))
	diag(yhat) = NA ; pgbmeProb = c(yhat) ; rm(yhat)
	ddesign$pgbmeProb = pgbmeProb[!is.na(pgbmeProb)]

	# add in fold ids
	diag(foldMat) = NA ; foldVec = c(foldMat)
	ddesign$fold = foldVec[!is.na(foldVec)]

	# generate glm out of samp preds
	glmBeta = coef( glm(
		formula(paste0('y ~ ', paste(names(ddesign)[3:11], collapse='+'))), 
		family='binomial', data=ddesign[ddesign$fold!=fold,]
		) )
	outDesignArray = cbind(1, ddesign[ddesign$fold==fold,names(ddesign)[3:11]])
	glmPreds = data.matrix(outDesignArray) %*% glmBeta

	# keep data only from fold left out
	ddesign = ddesign[ddesign$fold==fold,]	
	ddesign$glmProb = pnorm(glmPreds)
	return(ddesign)
})

#
preds = do.call('rbind', predList)
ddesign = preds

# 
predDfs = list(
	GLM = cbind(ddesign[,c('y','glmProb')], model='GLM'),
	GBME = cbind(ddesign[,c('y','gbmeProb')], model='GBME'),
	'P-GBME' = cbind(ddesign[,c('y','pgbmeProb')], model='P-GBME')
	)
predDfs = lapply(predDfs, function(x){names(x)[1:2]=c('actual','pred');return(x)})

# tabular data
aucSumm=do.call('rbind', lapply(predDfs,function(x){
	aucROC=roc.curve(x$pred[x$actual==1], x$pred[x$actual==0])$auc
	aucPR=pr.curve(x$pred[x$actual==1], x$pred[x$actual==0])$auc.integral
	return( c('AUC (ROC)'=aucROC,'AUC (PR)'=aucPR) ) }) )
aucSumm = aucSumm[order(aucSumm[,1],decreasing=TRUE),]
aucSumm = trim(format(round(aucSumm, 2), nsmall=2))
sink(paste0(path, 'table1_outSamplePerf.txt'))
print(aucSumm)
sink()
################################