#' sample from multivar norm
#' 
#' sample from multivar norm
#' @param mu mu
#' @param Sig2 Sig2

rmvnorm<-function(mu,Sig2){
  R<-t(chol(Sig2))
  return( R%*%(rnorm(length(mu),0,1)) + mu )
}