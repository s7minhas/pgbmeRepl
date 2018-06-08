# workspace ###############################
source('../setup.R')

# libs
libs = c('magic', 'msm', 'lme4', 'mnormt', 'abind', 'foreach', 'doParallel')
loadPkg(libs)
################################

# fns for p-gbme ###############################
pds <- 2012

# pull in mcmc code

source(paste0(funcPath, "gbme_partial.R"))

# load data
load(paste0(dataPath, 'modelData2012.rda'))

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
setwd(paste0(mainPath, 'inperf/'))
set.seed(6886)
est <- gbme.ef(
	y = y, Xd = xDyad, Xs = xNode, Xr = xNode, 
	k = 2, rho.calc = FALSE,
	NS = 3e+5, burn = 0, odens = 100, 
	out.name = "PGBME/out.2012"
	)
################################