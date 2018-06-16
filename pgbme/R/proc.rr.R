#' Procrustes transformation
#' 
#' Perform procrustes transformation
#' @param Y Y
#' @param X X

proc.rr<-function(Y,X){
  k<-dim(X)[2]
  A<-t(Y)%*%(X%*%t(X))%*%Y
  eA<-eigen(A,symmetric=T)
  Ahalf<-eA$vec[,1:k]%*%diag(sqrt(eA$val[1:k]),nrow=k)%*%t(eA$vec[,1:k])
  return( t(t(X)%*%Y%*%solve(Ahalf)%*%t(Y)) )
}