loadPkg(c('RColorBrewer'))

# Plot roc curves, depends RColorBrewer
# plot_type is "roc" or "pr"
rocPlot = function(rocData, type='roc', legPos=c(.56,.25),
  colorPal = 'Set1', colorManual=NULL, linetypes, legText=6, legSpace=3){

  if(type=='roc'){ 
    tmp=ggplot(rocData, aes(x=FPR, y=TPR, color=model, linetype=model)) + 
      geom_abline(intercept=0, slope=1, color='darkgrey') + 
      ylab('True Positive Rate (Sensitivity)') + xlab('False Positive Rate (1-Specificity)')

  }

  if(type=='pr'){ 
    tmp=ggplot(rocData, aes(x=rec, y=prec, color=model, linetype=model)) + 
      ylab('Precision') + xlab('Recall (True Positive Rate)')
  }

  if(is.null(colorManual)){
    tmp = tmp + scale_color_brewer(palette=colorPal)
  } else {
    tmp = tmp + scale_color_manual(values=colorManual)
  }

  tmp=tmp + 
    geom_line(lwd=1) +
    ylim(0,1) + 
    scale_linetype_manual(values=linetypes) + 
    theme(
      legend.position=legPos, legend.title=element_blank(),
      legend.background=element_blank(), 
      legend.text.align = 0, legend.text=element_text(size=legText),
      legend.key=element_rect(colour = NA, fill = NA), legend.key.size=unit(legSpace,'lines'),
      axis.ticks=element_blank(),    
      panel.border=element_blank()
    )
  return(tmp)
}

# gg separation plot
# thanks to http://www.peterhaschke.com/r/2013/04/22/SeparationPlot.html
ggSep = function(actual, proba, color, lty, actLineSize=2, fPath, save=TRUE){
  color = c('white',color)
  sepData = data.frame(actual, proba)
  sepData = sepData[order(sepData$proba),]
  tmp=ggplot(sepData) + 
    geom_rect(aes(xmin = 0, xmax = seq(length.out = length(actual)),
      ymin = 0, ymax = 1), fill = "transparent") +
    geom_linerange(aes(size=factor(actual), color = factor(actual),
      ymin = 0, ymax = 1, x = seq(length.out = length(actual))), alpha = 0.5) +
    geom_line(aes(y = proba, x = seq(length.out = length(actual)), linetype=lty), lwd = 4) + 
    scale_linetype_manual(values=lty) +
    scale_size_manual(values=c(1,actLineSize)) + 
    scale_color_manual(values=color) + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0), breaks=seq(0,1,.25)) + 
    theme(
      legend.position='none', 
      panel.grid=element_blank(), panel.border=element_rect(colour = "grey13"),
      axis.ticks=element_blank(),
      axis.text=element_blank(),
      axis.title=element_blank()
      )
  if(save){ ggsave(tmp, file=fPath, width=12, height=2) } else { return(tmp) }
}
####################################################################