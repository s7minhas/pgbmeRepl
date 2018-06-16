#' fisher.inv
#' 
#' fisher.inv
#' @param z z

fisher.inv <- function(z){
	return( (exp(2*z) - 1)/(exp(2*z) + 1) )
}