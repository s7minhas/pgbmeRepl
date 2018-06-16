#' Y.start
#' 
#' Y.start
#' @param y y

Y.start <- function(y){
  Y <- cbind(y, y)
  Y <- vect.mat(Y)
  diag(Y) <- 0
  return(Y)
}