#' covariance matrix for ab
#' 
#' covariance matrix for ab
#' @param a a 
#' @param b b
#' @param S0 S0
#' @param v0 v0

rSab.gibbs<-function(a,b,S0,v0){
  n<-length(a)
  ab<-cbind(a,b)
  Sn<-S0+ (t(ab)%*%ab)
  return( solve(rwish(solve(Sn),v0+n) ) )
}