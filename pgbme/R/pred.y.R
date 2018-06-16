#' Sample Y from posterior predictive distribution
#' 
#' Sample Y from posterior predictive distribution
#' @param theta theta
#' @param rho rho
#' @param se se
#' @param fam fam
#' @export

pred.y <- function(theta, rho, se, fam="gaussian"){
  n <- nrow(theta)
  e <- rmnorm(n*(n-1)/2, varcov = se*covmat(rho))
  e <- vect.mat(e)
  Y <- theta + e
  diag(Y) <- 0
  if (fam=="binomial") Y = 1*(Y>0)
  return(Y)
}