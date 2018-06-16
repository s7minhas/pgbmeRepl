#' fisher
#' 
#' fisher
#' @param rho rho

fisher <- function(rho){
	return(.5*log((1+rho)/(1-rho)))
}