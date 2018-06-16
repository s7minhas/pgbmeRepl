#' uv calc
#' 
#' uv calc
#' @param E E
#' @param UT UT

uv.E<-function(E, UT){
  u<- c(  t( (  E + t(E) )  *UT ) )
  u<-u[!is.na(u)]
  v<-c(  t( (  E - t(E) )  *UT ) )
  v<-v[!is.na(v)]
  return( list(u=u,v=v) )
}