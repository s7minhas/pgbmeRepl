# workspace ###############################
source('../setup.R')

#
loadPkg(c('reshape2', 'tidyr','dplyr'))
################################

# fns for gbme ###############################
# load data
load( paste0(dataPath, 'modelData2012.rda') )

# sample from posterior of imputed datasets
y <- data.matrix(bit.acc.t[[ '2012' ]])
cntryKey = data.frame(name=rownames(y), id=1:nrow(y), stringsAsFactors=FALSE)
yMelt = reshape2::melt(y) %>% dplyr::mutate(id=paste(Var1,Var2,sep='_'))

#
xDyad1 = xData[[10]]$xDyad
for(p in 1:dim(xDyad1)[3]){
	rownames(xDyad1[,,p])=colnames(xDyad1[,,p])=cntryKey$id
	diag(xDyad1[,,p])=NA }
xDyadMelt = reshape2::melt(xDyad1) %>% 
	tidyr::spread(key=Var3, value=value) %>%
	dplyr::mutate(id=paste(Var1,Var2,sep='_'))

#
xNode1 = data.frame(xData[[10]]$xNode,stringsAsFactors = FALSE)
xNode1$cntry = rownames(xNode1)
rm(xData)

# created stacked dyad matrix for glm
dyadData = yMelt[yMelt$Var1!=yMelt$Var2,]
# merge in dyad vars
dyadVars  = names(xDyadMelt)[3:7]
for(v in dyadVars){
	dyadData$tmp = xDyadMelt[match(dyadData$id,xDyadMelt$id),v]
	names(dyadData)[ncol(dyadData)] = v }
# merge in nodal vars
nodeVars = names(xNode1)[-ncol(xNode1)]
for(v in nodeVars){
	dyadData$tmp = xNode1[match(dyadData$Var1, xNode1$cntry),v]
	names(dyadData)[ncol(dyadData)] = paste0(v,'.sender')
	dyadData$tmp = xNode1[match(dyadData$Var2, xNode1$cntry),v]
	names(dyadData)[ncol(dyadData)] = paste0(v,'.receiver')	}

# mod spec
modForm = formula(
	paste0('value ~ ', 
		paste( names(dyadData)[5:ncol(dyadData)], collapse=' + ' )
		) )

# run glm
glmEst <- glm( modForm,
	data=dyadData, family=binomial(link='logit') )

# get preds
glmBeta = coef(glmEst)
glmPreds = data.matrix(cbind(1,dyadData[,names(glmBeta)[-1]])) %*% glmBeta
glmPreds = cbind(dyadData[,c('Var1','Var2')], glmProb = 1/(1+exp(-glmPreds)))

# save
setwd(paste0(mainPath,'inperf/GLM/'))
save(glmPreds, file='glmPreds.rda')
################################