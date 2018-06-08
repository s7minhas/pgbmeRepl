# workspace ###############################
source('../setup.R')

loadPkg(c('reshape2','ggplot2','magrittr', 'PRROC'))
theme_set(theme_bw())

# bin perf helpers
source(paste0(funcPath, 'binPerfHelpers.R'))
################################

################################
# load pred file
load(paste0(mainPath, 'inperf/modPerf.rda')) # ddesign
load(paste0(mainPath, 'inperf/GLM/glmPreds.rda')) # glmPreds

# merge glmPreds into ddesign
ddesign$id = paste(ddesign$Var1, ddesign$Var2, sep='_')
glmPreds$id = paste(glmPreds$Var1, glmPreds$Var2, sep='_')
ddesign$glmProb = glmPreds$glmProb[match(ddesign$id, glmPreds$id)]
################################

################################
# organize preds into list
predDfs = list(
	GLM = cbind(ddesign[,c('y','glmProb')], model='GLM'),
	GBME = cbind(ddesign[,c('y','gbmeProb')], model='GBME'),
	'P-GBME' = cbind(ddesign[,c('y','pgbmeProb')], model='P-GBME')
	)
predDfs = lapply(predDfs, function(x){names(x)[1:2]=c('actual','pred');return(x)})

# tabular data
aucSumm=do.call('rbind', lapply(predDfs,function(x){
	aucROC=roc.curve(x$pred[x$actual==1], x$pred[x$actual==0])$auc
	aucPR=pr.curve(x$pred[x$actual==1], x$pred[x$actual==0])$auc.integral
	return( c('AUC'=aucROC,'AUC (PR)'=aucPR) ) }) )
aucSumm = aucSumm[order(aucSumm[,1],decreasing=TRUE),]
aucSumm = trim(format(round(aucSumm, 2), nsmall=2))
print(aucSumm)
################################