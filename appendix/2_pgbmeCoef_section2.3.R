# workspace ###############################
rm(list=ls())
path <- '/home/minhas/' # ubuntu path format for ec2
# path <- '~/Research/pgbmeRepl/' # example path format for mac
aPath <- paste0(path, 'appendix/')
mPath <- paste0(path, 'main/')
setwd(aPath)

# install packages
toInstall <- c(
	'tidyr', 'reshape2', 'sbgcop', 'magic', 
	'msm', 'lme4', 'mnormt', 'abind', 'foreach', 
	'doParallel', 'ggplot2', 'dplyr', 
	'gridExtra', 'latex2exp')
for(pkg in toInstall){
  if(!pkg %in% installed.packages()[,1]){
    install.packages(pkg) } }

# load libraries
library(tidyr)
library(reshape2)
library(sbgcop)
library(magic)
library(msm)
library(lme4)
library(mnormt)
library(abind)
library(foreach)
library(doParallel)
library(ggplot2)
theme_set(theme_bw())
library(dplyr)
library(gridExtra)
library(latex2exp)

# helpers
char = function(x){as.character(x)}
num = function(x){as.numeric(char(x))}
cntr <- function(x) (x - mean(c(x), na.rm = TRUE))/sd(c(x), na.rm = TRUE)
trim = function (x) { gsub("^\\s+|\\s+$", "", x) }
################################

# load data files and assemble for imp ###############################
load(paste0(mPath, 'rawData.rda')) # bit.acc.t, covData, dist.norm

# only covData needs to be imputed
covDataL = do.call('rbind', covData)
forSbg = covDataL[,-c(1:3)]
sbgMod = sbgcop.mcmc(forSbg, nsamp=1000, odens=1, seed=6886)

# pull out a few imp datasets
set.seed(6886)
impIters = sample(501:1000, 11)
sbgImps = lapply(impIters, function(i){
	x = sbgMod$Y.impute[,,i] ; colnames(x) = colnames(sbgMod$Y.pmean)
	return( data.frame(covDataL[,1:3], x, row.names=NULL) ) })
sbgImps[[1]] = data.frame(covDataL[,1:3], sbgMod$Y.pmean,row.names=NULL)
################################

# design structs ###############################
pds <- char(1990:2012)
if(!all(paste0('results', pds, '.rda') %in% list.files())){
	cores <- length(pds)
	cl=makeCluster(cores) ; registerDoParallel(cl)
	shhh <- foreach(t = pds, 
		.packages = c('abind', 'magic', 'msm', 'lme4', 'mnormt', 'pgbme')
		) %dopar% {

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

		# organize data for pgbme
		y <- bit.acc.t[[ t ]]
		y <- apply(mat.vect(y), 1, prod)
		xDyadStart = xData[[1]]$xDyad
		xNodeStart = xData[[1]]$xNode
		impsToDraw = 2:length(xData)
		xDyad = lapply(xData[impsToDraw], function(x){x$xDyad})
		xNode = lapply(xData[impsToDraw], function(x){x$xNode})
		rm(xData)

	  	set.seed(6886)
		est <- pgbme(
			y = y, Xd = xDyadStart, Xs = xNodeStart, Xr = xNodeStart, 
			xInclImpList=TRUE, 
			Xd_L=xDyad, Xs_L=xNode, Xr_L=xNode,
			k = 2, rho.calc = FALSE,
			NS = 2e+4, burn = 1e+4, odens = 10
			)

		# save
		est <- est$est[,grepl('bd|bs|br',colnames(est$est))]
		colnames(est)[1:5] <- paste0('bd.',1:5)
		save(est, file=paste0('results',t,'.rda'))
	}
	stopCluster(cl)
}
################################

# var key ###############################
thetaL = lapply(pds, function(t){
	load(paste0('results',t,'.rda'))
	return( est ) })
varKey = data.frame( 
		var=c(
			paste0('bd.',1:5),
			paste0('bs',1:4),
			paste0('br',1:4)
			), stringsAsFactors = FALSE )
varKey$clean = c(
	'UDS Abs. Diff.$_{ij}$',
	'Law & Order Abs. Diff.$_{ij}$',
	'Log(GDP capita) Abs. Diff.$_{ij}$',
	'OECD$_{ij}$',
	'Distance$_{ij}$',
	'FDI/GDP_{i}', 'Disputes_{i}', 
	'GDP Capita Growth_{i}', 'PTAs_{i}',
	'FDI/GDP_{j}', 'Disputes_{j}',
	'GDP Capita Growth_{j}', 'PTAs_{j}'
	)

# helpful fns
facet_labeller = function(string){ TeX(string) }
summStats = function(x){ c(mu=mean(x), quantile(x, probs=c(0.025,0.05,0.95,0.975))) }
################################

# check convergence ###############################
ggConv = function(theta, xTitle, xLabAdd=FALSE){
	g = ggplot(theta, aes(x=imp, y=value, group=name)) + 
		geom_line() + 
		xlab(paste0('MCMC Chain for ', xTitle, ' Model')) + ylab('') + 
		facet_grid(name~type, scales='free_y', 
			labeller=as_labeller(facet_labeller, default=label_parsed)) + 
		theme(
			axis.ticks=element_blank(),
			panel.border=element_blank(),
			axis.text=element_text(size=7 ),
			axis.title.x=element_text(size=9, face='bold' ),			
			strip.text.x = element_text(size = 10, color='white' ),
			strip.text.y = element_text(size = 7, color='white' ),
			strip.background = element_rect(fill = "#525252", color='#525252')
			)
	if(!xLabAdd){g = g + xlab('')}
	return(g) }

# generate trace plot for 2012
theta=data.frame(thetaL[[length(pds)]])
theta$imp = 1:nrow(theta)
theta = melt(theta, id='imp')
theta$name = varKey$clean[match(theta$variable, varKey$var)]
theta$type = 'Dyadic Covariates'
theta$type[grepl('_{i}',theta$name,fixed=TRUE)] = 'Sender Covariates'
theta$type[grepl('_{j}',theta$name,fixed=TRUE)] = 'Receiver Covariates'	
theta$name = factor(theta$name, levels=varKey$clean[-1])

#
dCoef = ggConv(theta[grepl('ij',theta$name),], pds[length(pds)])
sCoef = ggConv(theta[grepl('_{i}',theta$name,fixed=T),], pds[length(pds)], TRUE)
rCoef = ggConv(theta[grepl('_{j}',theta$name,fixed=T),], pds[length(pds)])
dsrCoef = grid.arrange(dCoef, sCoef, rCoef, ncol=3)

ggsave(dsrCoef, width=10, height=7.5, file='figureA2.png')
################################

# parameter estimates over time ###############################
coefData = do.call('rbind', lapply(1:length(pds), function(i){
	theta=thetaL[[i]]
	summ = data.frame(t(apply(theta, 2, summStats)), stringsAsFactors=NULL)
	summ$var = rownames(summ) ; summ$name = varKey$clean[match(summ$var, varKey$var)]
	colnames(summ)[2:5] = c('lo95','lo90','up90','up95')
	summ$yr = pds[i] ; summ$sig = NULL
	summ$sig[summ$lo90 > 0 & summ$lo95 < 0] = "Positive at 90"
	summ$sig[summ$lo95 > 0] = "Positive"
	summ$sig[summ$up90 < 0 & summ$up95 > 0] = "Negative at 90"
	summ$sig[summ$up95 < 0] = "Negative"
	summ$sig[summ$lo90 < 0 & summ$up90 > 0] = "Insig"	
	summ = summ[summ$name!='Intercept',]
	summ$name = factor(summ$name, levels=varKey$clean[-1])
	return(summ)
}) )

#
coefData$type = 'Dyadic Covariates'
coefData$type[grepl('_{i}',coefData$name,fixed=TRUE)] = 'Sender Covariates'
coefData$type[grepl('_{j}',coefData$name,fixed=TRUE)] = 'Receiver Covariates'

# calc avg across time
coefData = coefData %>% dplyr::group_by(name) %>% dplyr::mutate(muT=mean(mu))

#
coefColors = c("Positive"=rgb(54, 144, 192, maxColorValue=255), 
	"Negative"= rgb(222, 45, 38, maxColorValue=255),
	"Positive at 90"=rgb(158, 202, 225, maxColorValue=255), 
	"Negative at 90"= rgb(252, 146, 114, maxColorValue=255),
	"Insig" = rgb(150, 150, 150, maxColorValue=255))
ggCoef = function(ggData){
	ggplot(ggData, aes(x=factor(yr),y=mu,color=sig)) +
		geom_hline(aes(yintercept=0), linetype=2, color = "black") +
		geom_hline(aes(yintercept=muT), size=2, color='grey70', alpha=.6) +	
		geom_point() + xlab('') + ylab('') + 
		geom_linerange(aes(ymin=lo90,ymax=up90),size=.7) +		
		geom_linerange(aes(ymin=lo95,ymax=up95),size=.3) +	
		scale_x_discrete('', breaks=seq(1990,2012,2)) +
		scale_color_manual(values=coefColors) + 
		facet_grid(name~type, scales='free_y', #nrow=3, ncol=3, 
			labeller=as_labeller(facet_labeller, default=label_parsed)) +
		theme(
			axis.ticks=element_blank(),
			panel.border=element_blank(),
			legend.position='none',
			axis.text.x=element_text(angle=45, hjust=1, size=7 ),
			axis.text.y=element_text(size=6 ),
			strip.text.x = element_text(size = 10, color='white' ),
			strip.text.y = element_text(size = 7, color='white' ),
			strip.background = element_rect(fill = "#525252", color='#525252')					
			)
	}

#
dCoef = ggCoef(coefData[grepl('ij',coefData$name),])
sCoef = ggCoef(coefData[grepl('_{i}',coefData$name,fixed=T),])
rCoef = ggCoef(coefData[grepl('_{j}',coefData$name,fixed=T),])
dsrCoef = grid.arrange(dCoef, sCoef, rCoef, ncol=3)

#
ggsave(dsrCoef, width=10, height=7.5, file='figureA3.pdf')
################################