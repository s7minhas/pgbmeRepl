#' Sample sender and receiver effects
#' 
#' Sample sender and receiver effects
#' @param u u
#' @param v v
#' @param su su
#' @param sv sv
#' @param piS.bd piS.bd
#' @param mu mu
#' @param Sab Sab
#' @param n n
#' @param XTu XTu
#' @param XTv XTv
#' @param tXTuXTu tXTuXTu
#' @param tXTvXTv tXTvXTv

rbeta.d.sr.gibbs <- function(u,v,su,sv,piS.bd,mu,Sab,n,XTu,XTv,tXTuXTu,tXTvXTv){
  del<-Sab[1,1]*Sab[2,2]-Sab[1,2]^2
  iSab<-rbind(cbind( diag(rep(1,n))*Sab[2,2]/del ,-diag(rep(1,n))*Sab[1,2]/del),
              cbind( -diag(rep(1,n))*Sab[1,2]/del,diag(rep(1,n))*Sab[1,1]/del) )
  rd<-dim(as.matrix(piS.bd))[1]
  
  if(dim(piS.bd)[1]>0){
    cov.beta.sr<-matrix(0,nrow=rd,ncol=2*n)
    iS<-rbind(cbind(solve(piS.bd),cov.beta.sr),
              cbind(t(cov.beta.sr),iSab)) }
  else{iS<-iSab}
  
  Sig <- chol2inv(chol(iS + tXTuXTu/su + tXTvXTv/sv))
  #this may have a closed form expression
  
  M<-Sig%*%(t((u%*%XTu)/su + (v%*%XTv)/sv) + iS%*%mu)
  
  beta.sr<-rmvnorm(M, Sig)
  return( list(beta.d=beta.sr[(rd>0):rd],s=beta.sr[rd+1:n],r=beta.sr[rd+n+1:n]) )
}