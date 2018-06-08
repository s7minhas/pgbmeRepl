# workspace ###############################
source('../setup.R')

# libs
libs = c('magic', 'msm', 'lme4', 'mnormt', 'abind', 'foreach', 'doParallel')
loadPkg(libs)
################################

# fns for p-gbme ###############################
pds <- 2012

# pull in mcmc code
source(paste0(funcPath, "gbme_partial_crossval.R"))

# load data
load( paste0(dataPath, 'modelData2012.rda') )

# sample from posterior of imputed datasets
y <- data.matrix(bit.acc.t[[ '2012' ]])
cntryKey = data.frame(cbind(name=rownames(y), id=1:nrow(y)))
rownames(y) = colnames(y) = cntryKey$id
y <- apply(mat.vect(y), 1, prod)

#
xDyad = xData[[10]]$xDyad
for(p in 1:dim(xDyad)[3]){ diag(xDyad[,,p])=NA }
dimnames(xDyad)[[1]] = dimnames(xDyad)[[2]] = cntryKey$id

#
xNode = xData[[10]]$xNode
rownames(xNode) = cntryKey$id
rm(xData)

# throw in some NAs to test ... set up k-val
set.seed(6886)
foldMat = matrix(sample(1:10, nrow(bit.acc.t[[ '2012' ]])^2, replace=TRUE), nrow=150)

# run folds in parallel
loadPkg(c('foreach','doParallel'))
cores = detectCores() - 1
cl=makeCluster(cores) ; registerDoParallel(cl)
foreach(fold = 1:10, .packages=c('lme4','magic','msm','mnormt')) %dopar% {

	#
	naMat = foldMat ; naMat[naMat==fold] <- NA ;
	naMat[naMat!=fold] <- 1 ; diag(naMat) = NA ;
	
	#
	yMiss = data.matrix(bit.acc.t[[ '2012' ]]) * naMat
	cntryKey = data.frame(cbind(name=rownames(yMiss), id=1:nrow(yMiss)))
	rownames(yMiss) = colnames(yMiss) = cntryKey$id
	yMiss <- apply(mat.vect(yMiss), 1, prod)

	#
	setwd(paste0(mainPath, 'outPerf/oPerf_PGBME/'))

	set.seed(6886)
	y=yMiss
	Xd = xDyad
	Xs = xNode
	Xr = xNode
	est <- gbme.ef(
		y = yMiss, Xd = xDyad, Xs = xNode, Xr = xNode, 
		k = 2, rho.calc = FALSE,
		NS = 1e+5, burn = 0, odens = 100, 
		out.name = paste("fold", fold, sep = "_")
		)
}
stopCluster(cl)
################################