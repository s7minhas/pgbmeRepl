# workspace ###############################
source('../setup.R')

#
loadPkg('lme4')
################################

# fns for gbme ###############################
pds <- 2012

# pull in mcmc code
source(paste0(mainPath, "functions/gbme.R"))

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

# throw in some NAs to test ... set up k-val
set.seed(6886)
foldMat = matrix(sample(1:10, nrow(y)^2, replace=TRUE), nrow=150)

# run folds in parallel
loadPkg(c('foreach','doParallel'))
cores = detectCores() - 1
cl=makeCluster(cores) ; registerDoParallel(cl)
foreach(fold = 1:10, .packages=c('lme4')) %dopar% {

	#
	naMat = foldMat ; naMat[naMat==fold] <- NA ;
	naMat[naMat!=fold] <- 1 ; diag(naMat) = NA ;
	yMiss = y * naMat

	#
	dir.create(paste0(mainPath, 'outPerf/oPerf_GBME/fold_',fold,'/'),showWarnings=FALSE)
	setwd(paste0(mainPath, 'outPerf/oPerf_GBME/fold_',fold,'/'))
	set.seed(6886)
	est <- gbme(
		Y = yMiss, Xd = xDyad1, Xs = xNode1, 
		fam='binomial', directed=FALSE, k = 2, 
		N = matrix(1,nrow(y),nrow(y)),
		NS = 1e+5, odens = 100, seed=6886,
		owrite=TRUE, awrite=TRUE
		)
}
stopCluster(cl)
################################