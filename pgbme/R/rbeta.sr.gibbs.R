#' Pull out multivariate dist of s/r
#' 
#' Pull out multivariate dist of s/r
#' @param s s
#' @param r r
#' @param X.u X.u
#' @param pim.b0sr pim.b0sr
#' @param piS.b0sr piS.b0sr
#' @param Sab Sab

rbeta.sr.gibbs<-function(s,r,X.u,pim.b0sr,piS.b0sr,Sab) {
  n <- length(s)
  del<-Sab[1,1]*Sab[2,2]-Sab[1,2]^2
  iSab<-rbind(cbind( diag(rep(1,n))*Sab[2,2]/del ,-diag(rep(1,n))*Sab[1,2]/del),
              cbind( -diag(rep(1,n))*Sab[1,2]/del,diag(rep(1,n))*Sab[1,1]/del) )
  
  S<-solve( solve(piS.b0sr) + t(X.u)%*%iSab%*%X.u )
  mu<-S%*%(  solve(piS.b0sr)%*%pim.b0sr+ t(X.u)%*%iSab%*%c(s,r))
  return(rmvnorm( mu,S))
}