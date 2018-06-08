# workspace ###############################
source('../setup.R')  

loadPkg(c('dplyr','reshape2','ggplot2','magrittr',
	'latex2exp','RColorBrewer'))
theme_set(theme_bw())
################################

################################
# load sim results
load(paste0(mainPath, 'simulation/monte_carlo_sims.Rdata'))
outBeta = lapply(out, function(x){
	est = x$'est'
	names(est)[2:3] = c("bd.1", "bd.2")
	return(est[c('Intercept',"bd.1", "bd.2", "bs1", "br1")]) })
outBeta = lapply(outBeta, function(x){do.call('cbind', x)})

# 
varKey = data.frame(
	dirty=c("bd.1", "bd.2", "bs1", "br1"),
	clean=c('$\\beta^{(d,1)}$', '$\\beta^{(d,2)}$', '$\\beta^{(s)}$', '$\\beta^{(r)}$'),
	actVal=c(1, -1/2, 0, 1/2),
	stringsAsFactors = FALSE)
################################

# sim bias analysis ###############################
meanParamVal = lapply(outBeta, 
	function(x){ apply(x, 2, mean) }) %>%
	do.call('rbind', .) %>% melt()

# get rid of Intercept
meanParamVal = meanParamVal[meanParamVal$Var2!='Intercept',]

# add clean labels
meanParamVal$label = varKey$clean[match(meanParamVal$Var2, varKey$dirty)]

# calc bias
meanParamVal$actVal = varKey$actVal[match(meanParamVal$Var2, varKey$dirty)]
meanParamVal$bias = meanParamVal$value - meanParamVal$actVal
with(meanParamVal, tapply(bias, Var2, mean))

# viz
biasPlot=ggplot(meanParamVal, aes(x=label, y=value, fill=label)) +
	geom_jitter(alpha=.5) +		
	geom_boxplot(outlier.alpha=.01,alpha=.7) +
	geom_hline(aes(yintercept=actVal,color=label)) +	
	ylab('Parameter Estimate') + xlab('') +
	scale_fill_discrete(
		'',
		labels=lapply(varKey$clean, TeX)
		) +
	guides(color=FALSE) +
	theme(
		panel.border=element_blank(),
		axis.ticks=element_blank(),
		axis.text.x=element_blank(),
		legend.position = 'top'
		)
ggsave(biasPlot, file=paste0(graphicsPath, 'simBias.pdf'), width=6, height=4)
################################

# sim coverage analysis ###############################
coverFN = function(x,y){ ( x[1]<y & y<x[2] )*1 }
coverParamVal = lapply(outBeta, function(x){
	bd1Qt = quantile(x[,'bd.1'], probs=c(0.025, 0.975))
	bd1Cover = coverFN(bd1Qt, varKey$actVal[varKey$dirty=='bd.1'])

	bd2Qt = quantile(x[,'bd.2'], probs=c(0.025, 0.975))
	bd2Cover = coverFN(bd2Qt, varKey$actVal[varKey$dirty=='bd.2'])

	bs1Qt = quantile(x[,'bs1'], probs=c(0.025, 0.975))
	bs1Cover = coverFN(bs1Qt, varKey$actVal[varKey$dirty=='bs1'])		

	br1Qt = quantile(x[,'br1'], probs=c(0.025, 0.975))
	br1Cover = coverFN(br1Qt, varKey$actVal[varKey$dirty=='br1'])			

	out=c('bd.1'=bd1Cover, 'bd.2'=bd2Cover, 'bs1'=bs1Cover, 'br1'=br1Cover)
	return(out)
}) %>% do.call('rbind', .) %>% data.frame()

# coverage stats
print(apply(coverParamVal, 2, mean))
################################