# workspace ###############################
source('../setup.R')

loadPkg(c('reshape2','ggplot2','magrittr'))
theme_set(theme_bw())

pasteVec = function(x,y){ as.vector( outer( x, y, paste0 ) ) }
################################

# fns for p-gbme ###############################
pds <- 2012

# pull in mcmc code
source(paste0(funcPath, "gbme.R"))

# load data
load( paste0(dataPath, 'modelData2012.rda') )

# sample from posterior of imputed datasets
y <- data.matrix(bit.acc.t[[ '2012' ]])
xDyadStart = xData[[1]]$xDyad
xNodeStart = xData[[1]]$xNode
rm(xData)

# read in model results
gbmePath = paste0(mainPath, 'inperf/GBME/')
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
load(paste0(mainPath, 'inperf/PGBME/out.2012/yhat.Rdata'))
yhat <- apply(pnorm(yhat), 1, mean)
yhat <- matrix(yhat, sqrt(length(yhat)), sqrt(length(yhat)))
diag(yhat) = NA ; pgbmeProb = c(yhat)
ddesign$pgbmeProb = pgbmeProb[!is.na(pgbmeProb)]

# 
save(ddesign, file=paste0(mainPath, 'inperf/modPerf.rda'))
################################