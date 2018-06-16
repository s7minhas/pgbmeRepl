#' sample from wishart
#' 
#' sample from wishart
#' @param S0 S0
#' @param nu nu

rwish<-function(S0,nu){ 
  S<-S0*0
  for(i in 1:nu){ z<-rmvnorm(rep(0,dim(as.matrix(S0))[1]), S0)
                  S<-S+z%*%t(z)  }
  return( S )
}