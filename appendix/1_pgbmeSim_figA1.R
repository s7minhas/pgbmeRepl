# workspace ###############################
rm(list=ls())
path <- '/home/minhas/' # ubuntu path format for ec2
# path <- '~/Research/pgbmeRepl/' # example path format for mac
aPath <- paste0(path, 'appendix/')
mPath <- paste0(path, 'main/')
setwd(aPath)

# install packages
toInstall <- c(
	'magic', 'msm', 'lme4', 'mnormt', 'abind', 
	'parallel', 'reshape2', 'ggplot2', 
	'latex2exp', 'magrittr')
for(pkg in toInstall){
  if(!pkg %in% installed.packages()[,1]){
    install.packages(pkg) } }

# load libraries
library(magic)
library(msm)
library(lme4)
library(mnormt)
library(abind)
library(parallel)
library(reshape2)
library(ggplot2)
theme_set(theme_bw())
library(latex2exp)
library(magrittr)

# helpers
library(pgbme)
char = function(x){as.character(x)}
num = function(x){as.numeric(char(x))}
cntr <- function(x) (x - mean(c(x), na.rm = TRUE))/sd(c(x), na.rm = TRUE)
trim = function (x) { gsub("^\\s+|\\s+$", "", x) }
################################

# PARTIAL PROBIT MODEL SIMULATION ###############################
if(!file.exists('pgbme_sim_results.rda')){
	monte_carlo <- function(w){

		library(pgbme)

		# NETWORK SIMULATION
		n  <<- 100
		X1 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
		X2 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
		z <- as.matrix(rnorm(n)) # nodal predictor
		abSigma <- function(rho){ return(matrix(c(1, rho, rho, 1), 2, 2)) }
		P <- rmnorm(n, varcov = abSigma(0)) # Sender-Receiver random effects
		b <<- c(1, -1/2, 0, 1/2)
		S <- matrix(P[,1] + b[3]*z, n, n, byrow = FALSE)
		R <- matrix(P[,2] + b[4]*z, n, n, byrow = TRUE)
		U <- rmnorm(n, varcov = diag(1))
		V <- rmnorm(n, varcov = diag(1))
		Z <- U%*%t(V)

		# Sample a network
		Y <- pred.y(-1.5 + b[1]*X1 + b[2]*X2 + S + R + Z, 
			rho = 0, se = 1, fam = "binomial")
		y <- apply(mat.vect(Y), 1, prod)

		Xd <- array(0, dim = c(n, n, 2))
		Xd[,,1] <- X1
		Xd[,,2] <- X2

		est <- pgbme(
			y = y, Xd = Xd, Xs = z, Xr = z, 
			k = 1, rho.calc = FALSE,
			NS = 2e+4, burn = 1e+4, odens = 10, 
			xInclImpList=FALSE)
		est <- est$est[,c(2:3,5:6)]
		colnames(est)[1:2] <- c('bd.1','bd.2')
		return(est)
	}

	# Calculate the number of cores
	cores <- 35
	# Initiate cluster
	cl <- makeCluster(cores)
	clusterExport(cl, c("mPath","monte_carlo"))
	out <- parLapply(cl, 1:100, monte_carlo)
	stopCluster(cl)
	save(out, file='pgbme_sim_results.rda') ; rm(out)
}
load('pgbme_sim_results.rda')
################################

################################
# var key
varKey = data.frame(
	dirty=c("bd.1", "bd.2", "bs1", "br1"),
	clean=c('$\\beta^{(d,1)}$', '$\\beta^{(d,2)}$', '$\\beta^{(s)}$', '$\\beta^{(r)}$'),
	actVal=c(1, -1/2, 0, 1/2),
	stringsAsFactors = FALSE)
################################

# sim bias analysis ###############################
meanParamVal = lapply(out, 
	function(x){ apply(x, 2, mean) }) %>%
	do.call('rbind', .) %>% melt()

# add clean labels
meanParamVal$label = varKey$clean[match(meanParamVal$Var2, varKey$dirty)]

# calc bias
meanParamVal$actVal = varKey$actVal[match(meanParamVal$Var2, varKey$dirty)]
meanParamVal$bias = meanParamVal$value - meanParamVal$actVal

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
ggsave(biasPlot, file='figureA1.pdf', width=6, height=4)
################################

# sim coverage analysis ###############################
coverFN = function(x,y){ ( x[1]<y & y<x[2] )*1 }
coverParamVal = lapply(out, function(x){
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
sink(file='simCoverStats_sectionA1.txt')
print(apply(coverParamVal, 2, mean))
sink()
################################