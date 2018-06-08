#### PARTIAL PROBIT MODEL SIMULATION
monte_carlo <- function(w){
  
setwd("~/Research/pgbmeRepl/main/")
source("pgbme.R")
# source("~/Desktop/gbme_partial.R")

# NETWORK SIMULATION
n  <<- 200
X1 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
X2 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
X3 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
X4 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
X5 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
z <- as.matrix(rnorm(n)) # nodal predictor
P <- rmnorm(n, varcov = covmat(0)) # Sender-Receiver random effects
b <<- c(1, -1/2, 0, 1/2)
S <- matrix(P[,1] + b[3]*z, n, n, byrow = FALSE)
R <- matrix(P[,2] + b[4]*z, n, n, byrow = TRUE)
U <- rmnorm(n, varcov = diag(2))
V <- rmnorm(n, varcov = diag(2))
Z <- U%*%t(V)

# Sample a network
Y <- pred.y(-1.5 + b[1]*X1 + b[2]*X2 + S + R + Z, rho = 0, se = 1, fam = "binomial")
y <- apply(mat.vect(Y), 1, prod)

Xd <- array(0, dim = c(n, n, 5))
Xd[,,1] <- X1
Xd[,,2] <- X2
Xd[,,3] <- X3
Xd[,,4] <- X4
Xd[,,5] <- X5

est <- pgbme(
	y = y, Xd = Xd, Xs = z, Xr = z, k = 2, 
	NS = 1000, burn = 1, odens = 1, rho.calc = FALSE
	)
est
}

Sys.time()
m <- monte_carlo(1)
Sys.time()

# Check size 
save(m, file = "m.Rdata")


y.hat <- calc_yhat(m)

# plot(c(y.hat), c(m$'yhat'))

#### PARTIAL PROBIT MODEL SIMULATION
monte_carlo <- function(w){
  
setwd("~/Research/pgbmeRepl/main/")
source("pgbme.R")
# source("~/Desktop/gbme_partial.R")

# NETWORK SIMULATION
n  <<- 200
X1 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
X2 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
X3 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
X4 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
X5 <- rmnorm(n, varcov = diag(n)) # dyad-level predictor
z <- as.matrix(rnorm(n)) # nodal predictor
P <- rmnorm(n, varcov = covmat(0)) # Sender-Receiver random effects
b <<- c(1, -1/2, 0, 1/2)
S <- matrix(P[,1] + b[3]*z, n, n, byrow = FALSE)
R <- matrix(P[,2] + b[4]*z, n, n, byrow = TRUE)
U <- rmnorm(n, varcov = diag(2))
V <- rmnorm(n, varcov = diag(2))
Z <- U%*%t(V)

# Sample a network
Y <- pred.y(-1.5 + b[1]*X1 + b[2]*X2 + S + R + Z, rho = 0, se = 1, fam = "binomial")
y <- apply(mat.vect(Y), 1, prod)

Xd <- array(0, dim = c(n, n, 5))
Xd[,,1] <- X1
Xd[,,2] <- X2
Xd[,,3] <- X3
Xd[,,4] <- X4
Xd[,,5] <- X5

Xd2 = Xd + array(rnorm(n*n*2), dim=c(n,n,5))
Xd3 = Xd + array(rnorm(n*n*2), dim=c(n,n,5))
xdl = list(Xd, Xd2, Xd3)
zl = list(z, z, z)


est <- pgbme(
	y = y, Xd = Xd, Xs = z, Xr = z, k = 2, 
	NS = 1000, burn = 1, odens = 1, rho.calc = FALSE,
	xInclImpList=TRUE, Xd_L=xdl, Xs_L=zl, Xr_L=zl
	)
est
}

Sys.time()
m <- monte_carlo(1)
Sys.time()

# Check size 
save(m, file = "m2.Rdata")


y.hat <- calc_yhat(m, xInclImpList=TRUE)

# plot(c(y.hat), c(m$'yhat'))