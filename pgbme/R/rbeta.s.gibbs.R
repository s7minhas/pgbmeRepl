#' Pull out multivariate dist of s/r
#' 
#' Pull out multivariate dist of s/r
#' @param s s
#' @param X.u X.u
#' @param pim.b0s pim.b0s
#' @param piS.b0s piS.b0s
#' @param s2a s2a

rbeta.s.gibbs<-function(s,X.u,pim.b0s,piS.b0s,s2a) {
  iSa<-diag(rep(1/s2a,n))
  S<-solve( solve(piS.b0s) + t(X.u)%*%iSa%*%X.u )
  mu<-S%*%(  solve(piS.b0s)%*%pim.b0s+ t(X.u)%*%iSa%*%s)
  return( rmvnorm( mu,S) )
}