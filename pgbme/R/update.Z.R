#' Update Z
#' 
#' Update Z
#' @param y y
#' @param mu mu
#' @param rho rho
#' @param Z Z

update.Z <- function(y, mu, rho, Z){
  tau <- c(-Inf, 0, Inf) 
  Z  <- mat.vect(Z)
  mu <- mat.vect(mu)
  sd <- sqrt(1 - rho^2)
 
for (t in sample(1:2, 2)){
    w  <- which(y == 1)
    Z[w, t] <- rtnorm(
      length(w), mu[w, t] + rho*(Z[w, 3-t] - mu[w, 3-t]), 
      sd = sd, lower = 0)
    Z[w, 3-t] <- rtnorm(
      length(w), mu[w, 3-t] + rho*(Z[w, t] - mu[w, t]), 
      sd = sd, lower = 0)
    
    w  <- which(y == 0)
    Z[w, t] <- rtnorm(
      length(w), mu[w, t] + rho*(Z[w, 3-t] - mu[w, 3-t]), 
      sd = sd, upper = tau[(Z[w, 3-t] < 0) + 2])
    Z[w, 3-t] <- rtnorm(
      length(w), mu[w, 3-t] + rho*(Z[w, t] - mu[w, t]), 
      sd = sd, upper = tau[(Z[w, t] < 0) + 2])
  }
  
  return(vect.mat(Z))
  
}