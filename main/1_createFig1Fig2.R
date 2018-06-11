# workspace ###############################
rm(list=ls())
path <- '~/Research/pgbmeRepl/main/'
setwd(path)

# load libraries
library(magic)
library(msm)
library(lme4)
library(mnormt)
library(abind)
library(foreach)
library(doParallel)
library(reshape2)
library(ggplot2)
theme_set(theme_bw())

# helpers
source("pgbme.R")
char = function(x){as.character(x)}
num = function(x){as.numeric(char(x))}
cntr <- function(x) (x - mean(c(x), na.rm = TRUE))/sd(c(x), na.rm = TRUE)
trim = function (x) { gsub("^\\s+|\\s+$", "", x) }
################################

# fns for p-gbme ###############################
pds <- c(1995,2010)
mcmc <- function(t, year=pds) {
	# pull in mcmc code
	source("pgbme.R")

	# load data
	load( paste0('modelData',year[t],'.rda') )

	# organize data for pgbme
	y <- bit.acc.t[[ char(year[t]) ]]
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
		NS = 2e+4, burn = 1e+4, odens = 10,
		seed=6886
		)
	
	# save results
	save(est, file=paste0('results',year[t],'.rda'))
	}
################################

# run p-gbme ###############################
if(
  !all(
    paste0('results',c(1995,2010),'.rda') %in% list.files()
    )
  ){
  results <- lapply(1:length(pds), function(i){
    mcmc(i) }) 
}

# pull out yhats
yhats <- lapply(pds, function(yr){
  load(paste0('results',yr,'.rda'))
  return(calc_yhat(est,TRUE)) })
names(yhats) <- pds
################################

# figure 1  ###############################
load('rawData.rda')
chn95grey = dplotGG(
	yhats$'1995', bit.acc.t$'1995',
	"CHN", actual =T, countryLabel='China') +
  scale_color_manual(values=c('#969696', '#000000'))
chn10grey = dplotGG(
	yhats$'2010', bit.acc.t$'2010',
	"CHN", actual =T, countryLabel='China') +
  scale_color_manual(values=c('#969696', '#000000'))
ggsave(chn95grey, height=4, width=4, file='figure1_a.pdf')
ggsave(chn10grey, height=4, width=4, file='figure1_b.pdf')
################################

# figure 2 ###############################
# org data for usa10
usa10 = predict.gbme(yhats$'2010', bit.acc.t$'2010',
	TRUE, threshold=NULL) ; usa10$y = data.matrix(usa10$y)
predDF = data.frame(senProb=usa10$yhat['USA',], recProb=usa10$yhat[,'USA'], 
  cntry=rownames(usa10$yhat), row.names=NULL, stringsAsFactors = FALSE)
trueDF = data.frame(actBIT=usa10$y['USA',], cntry=rownames(usa10$y))
predDF$actBIT = trueDF$actBIT[match(predDF$cntry,trueDF$cntry)]
predDF = predDF[predDF$cntry!='USA',]
predDF$prob = (predDF$senProb+predDF$recProb)/2

# subset to top 20% by us send
predDF = predDF[order(predDF$recProb,decreasing = TRUE),]
predDF = predDF[predDF$recProb>=quantile(predDF$recProb, .9) | predDF$recProb<=quantile(predDF$recProb, .1),]
predDF$cntry = factor(predDF$cntry, levels=rev(predDF$cntry))

#
predDF$cntryTop = char(predDF$cntry) ; predDF$cntryBot = char(predDF$cntry)
for(i in 1:nrow(predDF)){
  if(i %% 2 == 0){predDF$cntryTop[i]=''}else{predDF$cntryBot[i]=''} }

# 
predDF$actBIT[predDF$actBIT==1] = 'BIT signed by 2010'
predDF$actBIT[predDF$actBIT==0] = 'No BIT signed by 2010'
ggcols = c('BIT signed by 2010'='#4393c3','No BIT signed by 2010'='#d6604d')

#
ggVertUSA = ggplot(predDF, aes(x=cntry, color=factor(actBIT))) +
  geom_linerange(aes(ymin=0,ymax=-recProb), linetype='solid', size=1) +
  geom_linerange(aes(ymax=0,ymin=senProb), linetype='solid', size=1) +
  geom_hline(aes(yintercept=0), color='gray40', size=1) +    
  geom_point(aes(y=-recProb), shape=15, size=1.7) +
  geom_point(aes(y=senProb), shape=17, size=1.7) +  
  geom_vline(aes(xintercept=16.5), size=.8, linetype='dashed', color='gray60') +      
  scale_y_continuous(breaks=seq(-1,1,.25), labels=abs(seq(-1,1,.25))) +
  scale_color_manual(values=c('#000000', '#969696')) +
  xlab('') + 
  ylab('Probability USA is demanded    Probability USA demands   \n   as treaty partner (square)     treaty from other (triangle)') +
  coord_flip() +
  theme(
    axis.ticks=element_blank(),
    panel.border=element_blank(),
    axis.text.x=element_text(size=7),
    axis.text.y=element_text(size=6),
    legend.position='bottom',
    legend.title=element_blank()
    )
ggsave(ggVertUSA, file='figure2.pdf', height=6, width=5)
################################