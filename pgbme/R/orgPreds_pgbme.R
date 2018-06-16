#' Organize predictions from pgbme
#' 
#' Organize predictions from pgbme
#' @param yhat Vector of predicted probabilities
#' @param y Actual values
#' @param directional logical
#' @param threshold keep this set to NULL
#' @export

orgPreds_pgbme <- function(yhat, y, directional = TRUE, threshold = NULL){
  yhat <- apply(pnorm(yhat), 1, mean)
  yhat <- matrix(yhat, sqrt(length(yhat)), sqrt(length(yhat)))
  colnames(yhat) <- rownames(yhat) <- rownames(y)
  diag(yhat) <- 0

if (!is.null(threshold)){
if (threshold == "optimal"){
  yobs    <- apply(mat.vect(y), 1, prod)
  cat("Calculating optimal threshold", "\n")
  opthres <- function(t){
    yhat <- mat.vect(yhat) 
    yhat <- 1*(yhat > t)
    yhat <- apply(yhat, 1, prod)
    - (
      sensitivity(
        as.factor(yhat), as.factor(yobs)) + specificity(as.factor(yhat), as.factor(yobs))
      )
  }
}
  
  threshold <- optim(0.5, opthres, method = "L-BFGS-B", lower = 0.01, upper = 0.99)$par
  cat(round(threshold, 2), "\n")
}
  
if (directional){
  if (!is.null(threshold)) yhat <- 1*(yhat > threshold)
    return(list(yhat = yhat, y = y))
}

if (!directional){
  
  if (!is.null(threshold)){
    yhat <- mat.vect(yhat) 
    yhat <- 1*(yhat > threshold)
    yhat <- apply(yhat, 1, prod)
  }
  
  if (is.null(threshold)){
    yhat <- mat.vect(yhat)
    yhat <- apply(yhat, 1, prod)
    yhat <- vect.mat(cbind(yhat, yhat))
    colnames(yhat) <- rownames(yhat) <- rownames(y)
  }

  return(list(yhat = yhat, y = y))
  }
}