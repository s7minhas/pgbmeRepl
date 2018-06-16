#' Calculate a+b+u'v
#' 
#' Calculate a+b+u'v
#' @param a a
#' @param b b
#' @param E E
#' @param e e
#' @param f f
#' @export

reef<-function(a,b,E,e,f){
  n<-length(a)
  return( matrix(a,nrow=n,ncol=n,byrow=F)+matrix(b,nrow=n,ncol=n,byrow=T)+E+e%*%t(f) )
}