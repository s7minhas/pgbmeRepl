#' Update y
#' 
#' Update y
#' @param y y
#' @param Y Y
#' @param m m
#' @param rho rho
#' @param theta theta

update.Y <- function(y, Y, m, rho, theta){
  z <- mat.vect(theta)
  Y[y==1, ] <- 1
  
  sampY <- function(t){
    Y[y==0 & Y[, 3-t] == 1, t] <- 0
    w <- which(y==0 & Y[, 3-t] == 0)
    Y[w, t] <- rbinom(
      length(w), 1, 
      pnorm(m[w, t] + rho*(z[w, 3-t] - m[w, 3-t]), sd = sqrt(1-rho^2)))
    Y[  ,t]
  }
  
  for (t in sample(1:2, 2)) Y[,t] <- sampY(t)
  
  return(Y)
}