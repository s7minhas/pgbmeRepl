# workspace ###############################
source('setup.R')
setwd(mainPath)

# libs
libs = c('reshape2', 'ggplot2', 'tidyr', 'dplyr', 'gridExtra', 'latex2exp')
loadPkg(libs)
theme_set(theme_bw())
################################

# useful params ###############################
load( paste0(dataPath, 'modelData1990.rda') )
varKey = data.frame( cbind(
		c(
			'Intercept',
			dimnames(xData[[1]]$xDyad)[[3]],
			rep(colnames(xData[[1]]$xNode),2)
			),
		c(
			'Intercept', 
			paste0('bd.',1:dim(xData[[1]]$xDyad)[3]),
			paste0('bs',1:ncol(xData[[1]]$xNode)),
			paste0('br',1:ncol(xData[[1]]$xNode))
			) ), stringsAsFactors = FALSE )
varKey$clean = 'Intercept'
varKey$clean[2:nrow(varKey)] = c(
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
rm(xData)

# bring in theta data
yrs = 1990:2012
thetaL = lapply(yrs, function(t){
	f=paste0(mainPath, 'MCMC/out.',t, '/theta')
	return( read.table(f, header=TRUE) ) })

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
theta=thetaL[[length(yrs)]]
theta=theta[
	1000:nrow(theta), 
	names(theta)[grepl('bd|bs|br|Intercept', names(theta))]]
theta$imp = 1:nrow(theta)
theta = melt(theta, id='imp')
theta$name = varKey$clean[match(theta$variable, varKey$X2)]
theta = theta[theta$name!='Intercept',]
theta$type = 'Dyadic Covariates'
theta$type[grepl('_{i}',theta$name,fixed=TRUE)] = 'Sender Covariates'
theta$type[grepl('_{j}',theta$name,fixed=TRUE)] = 'Receiver Covariates'	
theta$name = factor(theta$name, levels=varKey$clean[-1])

#
dCoef = ggConv(theta[grepl('ij',theta$name),], yrs[length(yrs)])
sCoef = ggConv(theta[grepl('_{i}',theta$name,fixed=T),], yrs[length(yrs)], TRUE)
rCoef = ggConv(theta[grepl('_{j}',theta$name,fixed=T),], yrs[length(yrs)])
dsrCoef = grid.arrange(dCoef, sCoef, rCoef, ncol=3)

ggsave(dsrCoef, width=10, height=7.5, file=paste0(graphicsPath, 'conv2012.png') )
################################

# parameter estimates over time ###############################
coefData = do.call('rbind', lapply(1:length(yrs), function(i){
	theta=thetaL[[i]]
	theta=theta[
		1000:nrow(theta), # burn more?
		names(theta)[grepl('bd|bs|br|Intercept', names(theta))]]
	summ = data.frame(t(apply(theta, 2, summStats)), stringsAsFactors=NULL)
	summ$var = rownames(summ) ; summ$name = varKey$clean[match(summ$var, varKey$X2)]
	colnames(summ)[2:5] = c('lo95','lo90','up90','up95')
	summ$yr = yrs[i] ; summ$sig = NULL
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
ggsave(dsrCoef, width=10, height=7.5, 
	file=paste0(graphicsPath, 'coefSumm.pdf') )
################################