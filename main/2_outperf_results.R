# workspace ###############################
source('../setup.R')

loadPkg(c('reshape2','ggplot2','magrittr'))
theme_set(theme_bw())

pasteVec = function(x,y){ as.vector( outer( x, y, paste0 ) ) }
################################

# quick and dirty cross-val summary ###############################
pds <- 2012

# pull in mcmc code
source(paste0(funcPath, "gbme.R"))

# load data
load( paste0(dataPath, 'modelData2012.rda') )

# sample from posterior of imputed datasets
y <- data.matrix(bit.acc.t[[ '2012' ]])
xDyadStart = xData[[10]]$xDyad
xNodeStart = xData[[10]]$xNode
rm(xData)

predList = lapply(1:10, function(foldID){

	# read in gbme model results
	gbmePath = paste0(mainPath, 'outPerf/oPerf_GBME/fold_',foldID,'/')
	OUT = read.table(paste0(gbmePath, 'OUT'), header=TRUE)
	a = read.table(paste0(gbmePath, 'A'), header=TRUE)
	Z = read.table(paste0(gbmePath, 'z'), header=TRUE)

	#convert to an array
	nss<-dim(OUT)[1]
	nss <- nrow(a)
	n<-dim(Z)[1]/nss
	k<-dim(Z)[2]
	PZ<-array(dim=c(n,k,nss))
	for(i in 1:nss) { PZ[,,i]<-as.matrix(Z[ ((i-1)*n+1):(i*n) ,])  }
	PZ<-PZ[,,-(1:round(nss/2))]     #drop first half for burn in

	#find posterior mean of Z %*% t(Z)
	ZTZ<-matrix(0,n,n)
	for(i in 1:dim(PZ)[3] ) { ZTZ<-ZTZ+PZ[,,i]%*%t(PZ[,,i]) }
	ZTZ<-ZTZ/dim(PZ)[3] 

	#a configuration that approximates posterior mean of ZTZ
	tmp<-eigen(ZTZ)
	Z.pm<-tmp$vec[,1:k]%*%sqrt(diag(tmp$val[1:k]))

	# burn first half of chain and get means
	OUT = OUT[round(nrow(OUT)/2,0):nrow(OUT),c('scan','b0',paste0('bd',1:5),paste0('bs',1:4))]
	betaTab = t(apply(OUT[,-1], 2, function(x){
		c(mu=mean(x), med=median(x), quantile(x, c(.025,.05,.95,.975))) }))
	rownames(betaTab) = c('Intercept', dimnames(xDyadStart)[[3]], colnames(xNodeStart))
	beta = apply(OUT, 2, mean)[-1]
	a = a[round(nrow(a)/2,0):nrow(a),]
	aMu = apply(a, 2, mean)

	# 
	ddesign = dcast(reshape2::melt(xDyadStart), Var1 + Var2  ~ Var3)
	for(v in colnames(xNodeStart)){
		ddesign$tmp = xNodeStart[,v][match(ddesign$Var1, rownames(xNodeStart))]
		names(ddesign)[ncol(ddesign)] = v }
	ddesign = ddesign[ddesign$Var1 != ddesign$Var2, ]
	ddesign$y = melt(y)$value[match(
		paste0(ddesign$Var1,'_',ddesign$Var2),
		paste0(melt(y)$Var1,'_',melt(y)$Var2) )]
	ddesign$a = aMu[match(ddesign$Var1, rownames(y))]
	tmp = Z.pm %*% t(Z.pm) ; diag(tmp) = NA
	ddesign$ztz = tmp[!is.na(tmp)] ; rm(tmp)

	#
	X = cbind(1, ddesign[,c(dimnames(xDyadStart)[[3]], colnames(xNodeStart))])
	gbmePred = data.matrix(X) %*% beta + ddesign$a + ddesign$ztz
	ddesign$gbmeProb = c(pnorm(gbmePred))

	#
	load(paste0(mainPath, 'outPerf/oPerf_PGBME/fold_',foldID,'/yhat.Rdata'))
	yhat <- apply(pnorm(yhat), 1, mean)
	yhat <- matrix(yhat, sqrt(length(yhat)), sqrt(length(yhat)))
	diag(yhat) = NA ; pgbmeProb = c(yhat)
	ddesign$pgbmeProb = pgbmeProb[!is.na(pgbmeProb)]

	# add in fold ids
	set.seed(6886)
	foldMat = matrix(sample(1:10, nrow(bit.acc.t[[ '2012' ]])^2, replace=TRUE), nrow=150)
	diag(foldMat) = NA ; foldVec = c(foldMat)
	ddesign$fold = foldVec[!is.na(foldVec)]

	# generate glm out of samp preds
	glmBeta = coef( glm(
		formula(paste0('y ~ ', paste(names(ddesign)[3:11], collapse='+'))), 
		family='binomial', data=ddesign[ddesign$fold!=foldID,]
		) )
	outDesignArray = cbind(1, ddesign[ddesign$fold==foldID,names(ddesign)[3:11]])
	glmPreds = data.matrix(outDesignArray) %*% glmBeta
	glmProbs = 1/(1+exp(-glmPreds))

	# keep data only from fold left out
	ddesign = ddesign[ddesign$fold==foldID,]	
	ddesign$glmProb = glmProbs
	return(ddesign)
})

#
preds = do.call('rbind', predList)
ddesign = preds

# 
predDfs = list(
	GLM = cbind(ddesign[,c('y','glmProb')], model='GLM'),
	GBME = cbind(ddesign[,c('y','gbmeProb')], model='GBME'),
	'P-GBME' = cbind(ddesign[,c('y','pgbmeProb')], model='P-GBME')
	)
predDfs = lapply(predDfs, function(x){names(x)[1:2]=c('actual','pred');return(x)})

# tabular data
loadPkg('PRROC')
aucSumm=do.call('rbind', lapply(predDfs,function(x){
	aucROC=roc.curve(x$pred[x$actual==1], x$pred[x$actual==0])$auc
	aucPR=pr.curve(x$pred[x$actual==1], x$pred[x$actual==0])$auc.integral
	return( c('AUC'=aucROC,'AUC (PR)'=aucPR) ) }) )
aucSumm = aucSumm[order(aucSumm[,1],decreasing=TRUE),]
aucSumm = trim(format(round(aucSumm, 2), nsmall=2))
print(aucSumm)
################################