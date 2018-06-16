#' Matrix to vector
#' 
#' Back-transform N*(N-1)x2 vector into NxN matrix
#' @param Y a vector transformed from a matrix using mat.vect.
#' @export

vect.mat <- function(Y){
  n = (sqrt((8*nrow(Y)+1))+1)/2
  M = matrix(0, n, n)
  M[t(lower.tri(M))] = Y[,1]
  M = t(M)
  M[upper.tri(M)] = Y[,2] 
  return(M)
}