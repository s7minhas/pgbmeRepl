# workspace ###############################
rm(list=ls())
if(Sys.info()['user'] %in% c('janus829','s7m')){
	source('~/Dropbox/Research/BiProbPO/PartialNet/functions/setup.R')
	dir.create(paste0(mainPath, "PartialNet/MCMC"), showWarnings=FALSE)
	setwd(paste0(mainPath, "PartialNet"))
}

if(Sys.info()['user'] %in% c('minhas')){
	source('/home/minhas/PartialNet/functions/setup.R')
	dir.create(paste0(mainPath, "PartialNet/MCMC"), showWarnings=FALSE)
	setwd(paste0(mainPath, 'PartialNet/'))
}

# libs
libs = c('magic', 'msm', 'lme4', 'mnormt', 'abind', 'foreach', 'doParallel')
loadPkg(libs)
################################

# fns for p-gbme ###############################
pds <- 1990:2012
mcmc <- function(t, year=pds, returnOutput=FALSE) {
	# pull in mcmc code
	source(paste0(mainPath, "PartialNet/functions/gbme_partial_wImputation.R"))

	# load data
	load( paste0(dataPath, 'modelData',year[t],'.rda') )

	# sample from posterior of imputed datasets
	y <- bit.acc.t[[ char(year[t]) ]]
	y <- apply(mat.vect(y), 1, prod)
	xDyadStart = xData[[1]]$xDyad
	xNodeStart = xData[[1]]$xNode
	set.seed(6886) ; impsToDraw = sample(2:length(xData), 100)
	xDyad = lapply(xData[impsToDraw], function(x){x$xDyad})
	xNode = lapply(xData[impsToDraw], function(x){x$xNode})
	rm(xData)

	est <- gbme.ef_wImp(
		y = y, Xd = xDyad, Xs = xNode, Xr = xNode, 
		XdStart=xDyadStart, XsStart=xNodeStart, XrStart=xNodeStart, 
		k = 2, rho.calc = FALSE,
		NS = 3e+5, burn = 1.5e+5, odens = 100, 
		out.name = paste("MCMC/out", year[t], sep = ".")
		)
	if(returnOutput){ return(est) } else { return(print(paste0('Period ', year[t], ' done'))) }
}
################################

# run p-gbme ###############################
cl <- makeCluster(4) ; registerDoParallel(cl)
foreach(i = 1:length(pds), 
	.packages=c('abind', 'magic', 'msm', 'lme4', 'mnormt')) %dopar% mcmc(i)
stopCluster(cl)
################################