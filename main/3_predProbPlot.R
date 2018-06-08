# workspace ###############################
source('setup.R')
source(paste0(funcPath, "gbme_partial_wImputation.R"))
setwd(mainPath)

# libs
libs = c('reshape2', 'ggplot2', 'tidyr', 'latex2exp', 
	'igraph', 'countrycode', 'data.table', 
	'verification', 'coda', 'OptimalCutpoints', 'caret'
	)
loadPkg(libs)
theme_set(theme_bw())
# if Rgraphviz is not installed, use bioClite to install it
# source("https://bioconductor.org/biocLite.R") ; biocLite("Rgraphviz")
library("Rgraphviz")
################################

# Predicted network probabilities of ntw ties ###############################
predict.gbme <- function(year, directional = TRUE, threshold = NULL){
  load(paste0(mainPath, '/MCMC/out.',year, '/yhat.Rdata'))
  load( paste0(dataPath, 'modelData',year,'.rda') )
  y    <- bit.acc.t[[char(year)]]
  yhat <- apply(pnorm(yhat), 1, mean)
  yhat <- matrix(yhat, sqrt(length(yhat)), sqrt(length(yhat)))
  colnames(yhat) <- rownames(yhat) <- rownames(y)
  diag(yhat) <- 0

if (!is.null(threshold)){
if (threshold == "optimal"){
  yobs    <- apply(mat.vect(y), 1, prod)
  cat("Calculating optimal threshold", "\n")
  opthres <- function(t){
    yhat <- mat.vect(yhat) 
    yhat <- 1*(yhat > t)
    yhat <- apply(yhat, 1, prod)
    - (sensitivity(as.factor(yhat), as.factor(yobs)) + specificity(as.factor(yhat), as.factor(yobs)))
  }
}
  
  threshold <- optim(0.5, opthres, method = "L-BFGS-B", lower = 0.01, upper = 0.99)$par
  cat(round(threshold, 2), "\n")
}
  
if (directional){
  if (!is.null(threshold)) yhat <- 1*(yhat > threshold)
    return(list(yhat = yhat, y = y))
}

if (!directional){
  
  if (!is.null(threshold)){
    yhat <- mat.vect(yhat) 
    yhat <- 1*(yhat > threshold)
    yhat <- apply(yhat, 1, prod)
  }
  
  if (is.null(threshold)){
    yhat <- mat.vect(yhat)
    yhat <- apply(yhat, 1, prod)
    yhat <- vect.mat(cbind(yhat, yhat))
    colnames(yhat) <- rownames(yhat) <- rownames(y)
  }

  return(list(yhat = yhat, y = y))
  }
}
################################

# violin viz ###############################
# org data for usa10
usa10 = predict.gbme(2010, TRUE, threshold=NULL) ; usa10$y = data.matrix(usa10$y)
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
    legend.title=element_blank(),
    )

ggVertUSA = ggHorizUSA + 
  coord_flip() +
  theme(
    axis.text.x=element_text(size=7),
    axis.text.y=element_text(size=6)
    ) +
  scale_color_manual(values=c('#000000', '#969696'))
ggsave(ggVertUSA, file=paste0(graphicsPath, 'USA2010_gg_violinV_90_10_grey.pdf'), height=6, width=5)
################################

# Plotting fn ###############################
# Predicted probabilities directional network:
# $yhat predicted probabilities
# $y    optimal classification  
dplotGG <- function(year, country, actual=FALSE, threshold=NULL, countryLabel=country){
  yhat = predict.gbme(year, directional=TRUE, threshold=threshold)
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
################################

# plot ###############################
chn95grey = dplotGG(1995, "CHN", actual =T, countryLabel='China') +
  scale_color_manual(values=c('#969696', '#000000'))
chn10grey = dplotGG(2010, "CHN", actual =T, countryLabel='China') +
  scale_color_manual(values=c('#969696', '#000000'))
ggsave(chn95grey, height=4, width=4, file=paste0(graphicsPath, 'CHN1995_gg_grey.pdf'))
ggsave(chn10grey, height=4, width=4, file=paste0(graphicsPath, 'CHN2010_gg_grey.pdf'))
################################