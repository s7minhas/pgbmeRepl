# workspace ###############################
source('../setup.R')

# libs
libs = c('magic', 'msm', 'lme4', 'mnormt', 'abind', 'foreach', 'doParallel')
loadPkg(libs)
################################

# fns for p-gbme ###############################
pds <- 2012

# pull in mcmc code
source(paste0(mainPath, "functions/gbme_partial.R"))

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

#
loadPkg(c('foreach','doParallel'))
cores = detectCores() - 1
cl=makeCluster(cores) ; registerDoParallel(cl)
foreach(rank = 1:5, .packages=c('lme4','magic','msm','mnormt')) %dopar% {
	setwd(paste0(mainPath, 'var_by_k/'))
	set.seed(6886)
	est <- gbme.ef(
		y = y, Xd = xDyad, Xs = xNode, Xr = xNode, 
		k = rank, rho.calc = FALSE,
		NS = 2e+5, burn = 0, odens = 100, 
		out.name = paste("k", rank, sep = "_")
		)
}
stopCluster(cl)
################################