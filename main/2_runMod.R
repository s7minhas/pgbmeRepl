# workspace ###############################
setwd('~/Research/pgbmeRepl/main')
library(magic)
library(msm)
library(lme4)
library(mnormt)
library(abind)
library(foreach)
library(doParallel)

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

	# sample from posterior of imputed datasets
	y <- bit.acc.t[[ char(year[t]) ]]
	y <- apply(mat.vect(y), 1, prod)
	xDyadStart = xData[[1]]$xDyad
	xNodeStart = xData[[1]]$xNode
	impsToDraw = 2:6
	xDyad = lapply(xData[impsToDraw], function(x){x$xDyad})
	xNode = lapply(xData[impsToDraw], function(x){x$xNode})
	rm(xData)

	est <- pgbme(
		y = y, Xd = xDyadStart, Xs = xNodeStart, Xr = xNodeStart, 
		xInclImpList=TRUE, 
		Xd_L=xDyad, Xs_L=xNode, Xr_L=xNode,
		k = 2, rho.calc = FALSE,
		# NS = 3e+5, burn = 1.5e+5, odens = 100
		NS = 1e+5, burn = 5e+4, odens = 10,
		seed=6886
		)
	save(est, file=paste0('results',year[t],'.rda'))
	return(est) }
################################

# run p-gbme ###############################
cl <- makeCluster(2) ; registerDoParallel(cl)
results <- foreach(i = 1:length(pds), 
	.packages=c('abind', 'magic', 'msm', 'lme4', 'mnormt')) %dopar% mcmc(i)
stopCluster(cl)

# pull out yhats
yhats <- lapply(results, function(x){calc_yhat(x,TRUE)})
names(yhats) <- pds

# clean up workspace
rm(results)
################################

# start pred prob analysis ###############################
library(reshape2)
library(ggplot2)
theme_set(theme_bw())
library(tidyr)
library(latex2exp)
library(igraph)
library(countrycode)
library(data.table)
library(verification)
library(coda)
library(OptimalCutpoints)
library(caret)
################################

# figure 1  ###############################
# Predicted probabilities directional network:
# $yhat predicted probabilities
# $y    optimal classification  
dplotGG <- function(
	yPred, yAct, country, 
	actual=FALSE, threshold=NULL, 
	countryLabel=country){
  yhat = predict.gbme(yPred, yAct, directional=TRUE, threshold=threshold)
  pred = yhat$yhat
  predDF = data.frame(senProb=pred[country,], recProb=pred[,country],
    cntry=rownames(pred), row.names=NULL, stringsAsFactors = FALSE)
  if(!actual){
    gg = ggplot(predDF, aes(x=senProb, y=recProb))
  }
  if(actual){
    true = data.matrix(yhat$y)
    trueDF = data.frame(actBIT=true[country,], cntry=rownames(true),
      row.names=NULL, stringsAsFactors=FALSE)
    predDF$actBIT = trueDF$actBIT[match(predDF$cntry,trueDF$cntry)]
    gg = ggplot(predDF, aes(x=senProb, y=recProb, color=factor(actBIT))) +
      # scale_color_manual(values=c('1'='#67a9cf','0'='#ef8a62','NA'='white'))
    scale_color_manual(values=c('1'='#4393c3','0'='#d6604d','NA'='white'))
  }
  gg = gg + 
    geom_abline(color='grey40', linetype='dashed') +
    geom_text(aes(label=cntry), size=2) +
    labs(
      x=paste("Probability", countryLabel, "demands treaty from the other"),
      y=paste("Probability", countryLabel, "is demanded as treaty partner")
      ) +
    theme(
      axis.text = element_text(size=6),
      axis.title = element_text(size=8),
      axis.ticks=element_blank(),
      panel.border=element_blank(),
      legend.position='none'
      )
    return(gg)
}

# load y matrix
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
for(i in 1:nrow(predDF)){ if(i %% 2 == 0){predDF$cntryTop[i]=''}else{predDF$cntryBot[i]=''} }

# 
predDF$actBIT[predDF$actBIT==1] = 'BIT signed by 2010'
predDF$actBIT[predDF$actBIT==0] = 'No BIT signed by 2010'
ggcols = c('BIT signed by 2010'='#4393c3','No BIT signed by 2010'='#d6604d')

#
ggHorizUSA = ggplot(predDF, aes(x=cntry, color=factor(actBIT))) +
  geom_linerange(aes(ymin=0,ymax=-recProb), linetype='solid', size=1) +
  geom_linerange(aes(ymax=0,ymin=senProb), linetype='solid', size=1) +
  geom_hline(aes(yintercept=0), color='gray40', size=1) +    
  geom_point(aes(y=-recProb), shape=15, size=1.7) +
  geom_point(aes(y=senProb), shape=17, size=1.7) +  
  geom_vline(aes(xintercept=16.5), size=.8, linetype='dashed', color='gray60') +      
  scale_y_continuous(breaks=seq(-1,1,.25), labels=abs(seq(-1,1,.25))) +
  scale_color_manual(values=ggcols) +
  xlab('') + 
  ylab('Probability USA is demanded    Probability USA demands   \n   as treaty partner (square)     treaty from other (triangle)') +
  theme(
    axis.ticks=element_blank(),
    panel.border=element_blank(),
    axis.text.x=element_blank(),
    legend.position='bottom',
    legend.title=element_blank()
    )

ggVertUSA = ggHorizUSA + 
  coord_flip() +
  theme(
    axis.text.x=element_text(size=7),
    axis.text.y=element_text(size=6)
    ) +
  scale_color_manual(values=c('#000000', '#969696'))
ggsave(ggVertUSA, file='figure2.pdf', height=6, width=5)
################################