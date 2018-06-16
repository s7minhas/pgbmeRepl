#' predicted values for undirected case
#' 
#' predicted values for undirected case
#' @param beta.d beta.d
#' @param X.d X.d
#' @param s s
#' @param r r
#' @param E E
#' @param z z

theta.betaX.d.srE.z<-function(beta.d,X.d,s,r,E,z){
  m<-dim(X.d)[3]
  mu<-matrix(0,nrow=length(s),ncol=length(s))
  if(m>0){for(l in 1:m){ mu<-mu+beta.d[l]*X.d[,,l] }}
  tmp<-mu+re(s,r,E,z)
  diag(tmp)<-0
  return(tmp)
}