#' rse.beta.d.gibbs
#' 
#' rse.beta.d.gibbs
#' @param g0 g0
#' @param g1 g1
#' @param x x
#' @param XTx XTx
#' @param s s
#' @param r r
#' @param beta.d beta.d

rse.beta.d.gibbs<-function(g0,g1,x,XTx,s,r,beta.d){
  n<-length(s)
  return(1/rgamma(1, g0+choose(n,2)/2,g1+.5*sum( (x-XTx%*%c(beta.d,s,r))^2 ) ))
}