#' covmat
#' 
#' covmat
#' @param rho rho

covmat <- function(rho){
	return(matrix(c(1, rho, rho, 1), 2, 2))
}