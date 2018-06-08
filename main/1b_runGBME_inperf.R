# workspace ###############################
source('../setup.R')

#
loadPkg('lme4')
################################

# fns for gbme ###############################
pds <- 2012

# pull in mcmc code
source(paste0(funcPath, "gbme.R"))

# load data
load( paste0(dataPath, 'modelData2012.rda') )

# sample from posterior of imputed datasets
y <- data.matrix(bit.acc.t[[ '2012' ]])
cntryKey = data.frame(cbind(name=rownames(y), id=1:nrow(y)))
rownames(y) = colnames(y) = cntryKey$id


#
xDyad1 = xData[[10]]$xDyad
for(p in 1:dim(xDyad1)[3]){
	rownames(xDyad1[,,p])=colnames(xDyad1[,,p])=cntryKey$id
	diag(xDyad1[,,p])=NA }

#
xNode1 = xData[[10]]$xNode ; rownames(xNode1) = cntryKey$id
rm(xData)

#
setwd(paste0(mainPath, 'inperf/GBME/'))
set.seed(6886)
est <- gbme(
	Y = y, Xd = xDyad1, Xs = xNode1, 
	fam='binomial', directed=FALSE, k = 2, 
	N = matrix(1,nrow(y),nrow(y)),
	NS = 3e+5, odens = 100, seed=6886,
	owrite=TRUE, awrite=TRUE
	)
################################