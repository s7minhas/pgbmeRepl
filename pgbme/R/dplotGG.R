#' dplotGG
#' 
#' dplotGG
#' @param yPred yPred
#' @param yAct yAct
#' @param country country
#' @param actual actual
#' @param threshold threshold
#' @param countryLabel countryLabel
#' @export

dplotGG <- function(
  yPred, yAct, country, 
  actual=FALSE, threshold=NULL, 
  countryLabel=country){
  yhat = orgPreds_pgbme(yPred, yAct, directional=TRUE, threshold=threshold)
  pred = yhat$yhat
  # predDF = data.frame(recProb=pred[country,], senProb=pred[,country], ### pay attention here
  predDF = data.frame(senProb=pred[country,], recProb=pred[,country], ### pay attention here    
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