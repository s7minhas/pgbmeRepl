#' Sample dyadic beta and nodal effects
#' 
#' Sample dyadic beta and nodal effects
#' @param u u
#' @param su su
#' @param piS.bd piS.bd
#' @param mu mu
#' @param s2a s2a
#' @param rd rd
#' @param XTu XTu
#' @param tXTuXTu tXTuXTu

rbeta.d.s.gibbs <- function(u,su,piS.bd,mu,s2a,rd,XTu,tXTuXTu){
  iSa<-diag(rep(1/s2a,n))
  
  if(dim(piS.bd)[1]>0){
    cov.beta.s<-matrix(0,nrow=rd,ncol=n)
    iS<-rbind(cbind(solve(piS.bd),cov.beta.s),
              cbind(t(cov.beta.s),iSa)) }
  else{iS<-iSa}
  
  Sig<-chol2inv(chol(iS + tXTuXTu/su))
  #this may have a closed form expression
  
  theta <- Sig%*%(t((u%*%XTu)/su) + iS%*%mu)
  
  beta.s <- rmvnorm(theta,Sig)
  return( list(beta.d=beta.s[(rd>0):rd],s=beta.s[rd+1:n]) )
}