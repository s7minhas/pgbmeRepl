#' TuTv
#' 
#' TuTv
#' @param n n

TuTv<-function(n){
  Xu<-Xv<-NULL
  for(i in 1:(n-1)){
    tmp<-tmp<-NULL
    if( i >1 ){ for(j in 1:(i-1)){ tmp<-cbind(tmp,rep(0,n-i)) } }
    tmp<-cbind(tmp,rep(1,n-i)) 
    tmpu<-cbind(tmp,diag(1,n-i)) ; tmpv<-cbind(tmp,-diag(1,n-i))
    Xu<-rbind(Xu,tmpu) ; Xv<-rbind(Xv,tmpv)
  }
  
  return(list(Tu=cbind(Xu,Xu),Tv=cbind(Xv,-Xv)))
}