# workspace ###############################
source('../setup.R')  

# libs
libs = c('foreach', 'doParallel', 'PRROC', 'magrittr')
loadPkg(libs)
################################

################################
yhats = lapply(1:5, function(k){
	load(paste0(mainPath, 'var_by_k/k_',k,'/yhat.Rdata'))
	return(yhat) })

loadPkg(c('foreach','doParallel'))
cores=5
cl = makeCluster(cores)
registerDoParallel(cl)
preds <- foreach(yhat = yhats) %dopar% {
	yhat <- apply(pnorm(yhat), 1, mean)
	yhat <- matrix(yhat, sqrt(length(yhat)), sqrt(length(yhat)))
	diag(yhat) = NA
	return(c(yhat))
} 
predDF = data.frame(do.call('cbind', preds))
names(predDF) = paste0('k_',1:5)

# load data
load( paste0(dataPath, 'modelData2012.rda') )
predDF$y = c(data.matrix(bit.acc.t[[ '2012' ]]))
predDF = na.omit(predDF)

# run perf checks
perfStats = lapply(names(predDF)[1:5], function(k){
	aucROC=roc.curve(predDF[predDF$y==1,k], predDF[predDF$y==0,k])$auc
	aucPR=pr.curve(predDF[predDF$y==1,k], predDF[predDF$y==0,k])$auc.integral
	return( c('AUC (ROC)'=aucROC,'AUC (PR)'=aucPR) )
}) %>% do.call('rbind', .) 
rownames(perfStats) = names(predDF)[1:5]
print(round(perfStats, 2))
################################