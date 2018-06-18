# workspace ###############################
rm(list=ls())
path <- '/home/minhas/main/' # ubuntu path format for ec2
# path <- '~/Research/pgbmeRepl/replication/main/' # example path format for mac
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

# helpers
library(pgbme)
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
yPGBME <- apply(mat.vect(y), 1, prod)

#
xDyad = xData[[1]]$xDyad
for(p in 1:dim(xDyad)[3]){ diag(xDyad[,,p])=NA }
dimnames(xDyad)[[1]] = dimnames(xDyad)[[2]] = cntryKey$id

#
xNode = xData[[1]]$xNode
rownames(xNode) = cntryKey$id
rm(xData)
################################

# run pgbme ###############################
set.seed(6886)
if(!'resultsPGBME_inPerf.rda' %in% list.files()){
	pgbmeEst <- pgbme(
		y = yPGBME, Xd = xDyad, Xs = xNode, Xr = xNode, 
		k = 2, rho.calc = FALSE,
		NS = 2e+4, burn = 1e+4, odens = 10,
		xInclImpList=FALSE, seed=6886
		)
	save(pgbmeEst, file=paste0(path,'resultsPGBME_inPerf.rda'))
}
load(paste0(path,'resultsPGBME_inPerf.rda'))
yhat <- calc_yhat(pgbmeEst, FALSE); rm(pgbmeEst)
yhat <- apply(pnorm(yhat), 1, mean)
yhat <- matrix(yhat, sqrt(length(yhat)), sqrt(length(yhat)))
diag(yhat) = NA ; pgbmeProb = c(yhat) ; rm(yhat)
################################

# run gbme ###############################
source('gbme.R')
set.seed(6886)
setwd(paste0(path, 'gbme_inPerf/'))
if(!'OUT' %in% list.files()){
	gbme(
		Y = y, Xd = xDyad, Xs = xNode, 
		fam='binomial', directed=FALSE, k = 2, 
		N = matrix(1,nrow(y),nrow(y)),
		NS = 2e+4, odens = 10, seed=6886,	
		owrite=TRUE, awrite=TRUE
		)	
}
ddesign <- calc_yhat_gbme(y, xDyad, xNode, 
	paste0(path, 'gbme_inPerf/'))
ddesign$pgbmeProb = pgbmeProb[!is.na(pgbmeProb)]
################################

# run glm ###############################
xDyadMat <- apply(xDyad, 3, c)
n <- nrow(y) ; s <- rep(1:n, n) ; r <- rep(1:n, each=n)
glmEst  <- glm(
	c(y) ~ xDyadMat + xNode[s,] + xNode[r,], 
	family=binomial(link=probit))

# get preds
glmPreds = cbind(1,xDyadMat,xNode[s,],xNode[r,]) %*% coef(glmEst)
ddesign$glmProb = pnorm(glmPreds[!is.na(glmPreds)])
################################

# calc auc stats ###############################
# organize preds into list
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
sink(paste0(path, 'table1_inSamplePerf.txt'))
print(aucSumm)
sink()
################################