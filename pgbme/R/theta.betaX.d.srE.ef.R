#' predicted values for directed case
#' 
#' predicted values for directed case
#' @param beta.d beta.d
#' @param X.d X.d
#' @param s s
#' @param r r
#' @param E E
#' @param e e
#' @param f f

theta.betaX.d.srE.ef<-function(beta.d,X.d,s,r,E,e,f){
  m<-dim(X.d)[3]
  mu<-matrix(0,nrow=length(s),ncol=length(s))
  if(m>0){for(l in 1:m){ mu<-mu+beta.d[l]*X.d[,,l] }}
  tmp<-mu+reef(s,r,E,e,f)
  diag(tmp)<-0
  return(tmp)
}