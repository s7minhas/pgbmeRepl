#' Calculate a+b+u'u
#' 
#' Calculate a+b+u'u
#' @param a a
#' @param b b
#' @param E E
#' @param z z
#' @export

re<-function(a,b,E,z){
  n<-length(a)
  return(matrix(a,nrow=n,ncol=n,byrow=F)+matrix(b,nrow=n,ncol=n,byrow=T)+E+z%*%t(z) )
}