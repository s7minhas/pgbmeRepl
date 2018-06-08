#### PARTIAL PROBIT MODEL SIMULATION
monte_carlo <- function(w){

	source('../setup.R')  
	setwd(paste0(mainPath, "simulation/"))
	source("gbme_partial.R")

	# NETWORK SIMULATION
	n  <<- 100
	X1 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
	X2 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
	z <- as.matrix(rnorm(n)) # nodal predictor
	P <- rmnorm(n, varcov = covmat(0)) # Sender-Receiver random effects
	b <<- c(1, -1/2, 0, 1/2)
	S <- matrix(P[,1] + b[3]*z, n, n, byrow = FALSE)
	R <- matrix(P[,2] + b[4]*z, n, n, byrow = TRUE)
	U <- rmnorm(n, varcov = diag(1))
	V <- rmnorm(n, varcov = diag(1))
	Z <- U%*%t(V)

	# Sample a network
	Y <- pred.y(-1.5 + b[1]*X1 + b[2]*X2 + S + R + Z, 
		rho = 0, se = 1, fam = "binomial")
	y <- apply(mat.vect(Y), 1, prod)

	Xd <- array(0, dim = c(n, n, 2))
	Xd[,,1] <- X1
	Xd[,,2] <- X2

	est <- gbme.ef(
		y = y, Xd = Xd, Xs = z, Xr = z, k = 1, 
		NS = 20000, burn = 10000, odens = 1, 
		out.name = "out", rho.calc = FALSE)

	return(list(est = est, Y = Y))
}

loadPkg('parallel')
# Calculate the number of cores
cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(cores)
clusterExport(cl, c("monte_carlo"))
out <- parLapply(cl, 1:100, monte_carlo)
save(out, file = paste0(mainPath, "simulation/monte_carlo_sims.Rdata"))
stopCluster(cl)