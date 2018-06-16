#' Predictions from pgbme
#' 
#' Takes output from pgbme and returns matrix of predicted probabilities
#' @param m output from pgbme
#' @param xInclImpList logical indicating whether imputed covariate data was used in estimation of pgbme object
#' @export

calc_yhat <- function(m, xInclImpList=FALSE){
  # Dyadic coefficients:
  b <- m$est[,grep("bd", colnames(m$est))]
  # Empty zero matrix
  E <- matrix(0, nrow = nrow(m$Xd), ncol = nrow(m$Xd))  
  # Calculate predictions and collapse
  if(!xInclImpList){
    y_calc <- sapply(1:nrow(b), function(i){
      c(theta.betaX.d.srE.ef(
        b[i, ], m$Xd, m$s[i, ], m$r[i, ], E*0, 
        m$e[, , i], m$f[, , i]
        )) })
  } else {
    y_calc <- sapply(1:nrow(b), function(i){
      c(theta.betaX.d.srE.ef(
        b[i, ], m$Xd_L[[ m$xdId[i] ]], m$s[i, ], m$r[i, ], E*0, 
        m$e[, , i], m$f[, , i]
        )) })    
  }
  return(y_calc)
}