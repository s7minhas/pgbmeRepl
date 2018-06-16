#' XuXv
#' 
#' XuXv
#' @param X X

XuXv<-function(X){
  n <- nrow(X)
  Xu<-Xv<-NULL
  if(dim(X)[3]>0){
    for(r in 1:dim(X)[3]){
      xu<-xv<-NULL
      for(i in 1:(n-1)){
        for(j in (i+1):n){ xu<-c(xu,X[i,j,r]+X[j,i,r])
                           xv<-c(xv,X[i,j,r]-X[j,i,r]) }}
      Xu<-cbind(Xu,xu)
      Xv<-cbind(Xv,xv)  } 
  }
  return(list(Xu=Xu,Xv=Xv))
}