# workspace ###############################
source('setup.R')

# libs
libs = c('tidyr','reshape2','sbgcop','ggplot2','doParallel','foreach') ; loadPkg(libs)
################################

# load data files and assemble for imp ###############################
load(paste0(dataPath, 'rawData.rda')) # bit.acc.t, covData, dist.norm

# only covData needs to be imputed
covDataL = do.call('rbind', covData)
forSbg = covDataL[,-c(1:3)]
sbgMod = sbgcop.mcmc(forSbg, nsamp=1000, odens=1, seed=6886)

# convergence check
corCheck = sbgMod$C.psamp
for(i in 1:dim(corCheck)[3]){ diag(corCheck[,,i])=NA ; sbgMod$C.psamp[,,i][upper.tri(corCheck[,,i])] = NA  }
corMatChain = melt(corCheck) ; corMatChain = na.omit(corMatChain)
corMatChain$id = paste(corMatChain$Var1, corMatChain$Var2, sep='_')
sbgConv = ggplot(corMatChain, aes(x=Var3, y=value)) +
	geom_line() + ylab('') + xlab('') + 
	facet_wrap(~id, scales='free_y') + 
	theme(
		axis.ticks=element_blank(), 
		axis.text = element_text(size=5),
		strip.text = element_text(size=6)
		)
ggsave(paste0(graphicsPath, 'sbgConv.pdf'), sbgConv, device='pdf', width=10, height=10)

# pull out imp datasets (burn first 500)
sbgImps = lapply(500:dim(sbgMod$Y.impute)[3], function(i){
	x = sbgMod$Y.impute[,,i] ; colnames(x) = colnames(sbgMod$Y.pmean)
	return( data.frame(covDataL[,1:3], x, row.names=NULL) ) })
sbgImps[[1]] = data.frame(covDataL[,1:3], sbgMod$Y.pmean,row.names=NULL)
################################

# design structs ###############################
cores = detectCores() - 1
cl=makeCluster(cores) ; registerDoParallel(cl)
foreach(t = names(bit.acc.t)) %dopar% {

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
		file=paste0(dataPath, 'modelData',t,'.rda'))
}
stopCluster(cl)
################################
