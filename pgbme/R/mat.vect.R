#' Convert a matrix into a vector
#' 
#' Convert a matrix into a vector
#' @param Y Y is an n x n matrix
#' @export

mat.vect <- function(Y){ ## transform NxN matrix into N*(N-1)x2 vector
  Y.l = t(Y)[upper.tri(t(Y))]
  Y.u = Y[upper.tri(Y)]
  Y.new=cbind(Y.l, Y.u)
  return(Y.new)
}