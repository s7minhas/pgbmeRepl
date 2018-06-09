# workspace ###############################
setwd('~/Research/pgbmeRepl/main')
library(magic)
library(msm)
library(lme4)
library(mnormt)
library(abind)

# helpers
source("pgbme.R")
char = function(x){as.character(x)}
num = function(x){as.numeric(char(x))}
cntr <- function(x) (x - mean(c(x), na.rm = TRUE))/sd(c(x), na.rm = TRUE)
trim = function (x) { gsub("^\\s+|\\s+$", "", x) }
################################

# get data in order ###############################
# load data
load(paste0(dataPath, 'modelData2012.rda'))

# 
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
est <- pgbme(
	y = yPGBME, Xd = xDyad, Xs = xNode, Xr = xNode, 
	k = 2, rho.calc = FALSE,
	# NS = 3e+5, burn = 0, odens = 100, 
	NS = 1e4, burn = 5e3, odens = 10, 	
	xInclImpList=FALSE, seed=6886
	)
################################

# run gbme ###############################
source('gbme.R')
setwd(paste0(mainPath, 'inperf/GBME/'))
set.seed(6886)
est <- gbme(
	Y = y, Xd = xDyad, Xs = xNode, 
	fam='binomial', directed=FALSE, k = 2, 
	N = matrix(1,nrow(y),nrow(y)),
	# NS = 3e+5, odens = 100, 
	NS = 1e4, burn = 5e3, odens = 10, 		
	owrite=TRUE, awrite=TRUE, seed=6886
	)
################################

# run glm ###############################
library(reshape2)
library(dplyr)
library(tidyr)

#
yMelt = reshape2::melt(y) %>% dplyr::mutate(id=paste(Var1,Var2,sep='_'))
xDyadMelt = reshape2::melt(xDyad) %>% 
	tidyr::spread(key=Var3, value=value) %>%
	dplyr::mutate(id=paste(Var1,Var2,sep='_'))
xNode = data.frame(xNode,stringsAsFactors = FALSE)	

# created stacked dyad matrix for glm
dyadData = yMelt[yMelt$Var1!=yMelt$Var2,]
# merge in dyad vars
dyadVars  = names(xDyadMelt)[3:7]
for(v in dyadVars){
	dyadData$tmp = xDyadMelt[match(dyadData$id,xDyadMelt$id),v]
	names(dyadData)[ncol(dyadData)] = v }
# merge in nodal vars
nodeVars = names(xNode)[-ncol(xNode)]
for(v in nodeVars){
	dyadData$tmp = xNode[match(dyadData$Var1, xNode$cntry),v]
	names(dyadData)[ncol(dyadData)] = paste0(v,'.sender')
	dyadData$tmp = xNode[match(dyadData$Var2, xNode$cntry),v]
	names(dyadData)[ncol(dyadData)] = paste0(v,'.receiver')	}

# mod spec
modForm = formula(
	paste0('value ~ ', 
		paste( names(dyadData)[5:ncol(dyadData)], collapse=' + ' )
		) )

# run glm
glmEst <- glm( modForm,
	data=dyadData, family=binomial(link='logit') )

# get preds
glmBeta = coef(glmEst)
glmPreds = data.matrix(cbind(1,dyadData[,names(glmBeta)[-1]])) %*% glmBeta
glmPreds = cbind(dyadData[,c('Var1','Var2')], glmProb = 1/(1+exp(-glmPreds)))
################################