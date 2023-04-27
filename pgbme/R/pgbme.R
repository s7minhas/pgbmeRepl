#' Estimate a PGBME model
#' 
#' This function provides a MCMC routine to estimate a PGBME model. The model is 
#' detailed in "Modeling Asymmetric Relationships from Symmetric Networks" 
#' (Rozenas, Minhas & Ahlquist 2018). 
#' 
#' @param y a symmetric and unweighted square matrix (n*n).
#' @param Xd a n*n*rd array representing dyad-specific predictors.
#' @param Xs a n*rs matrix of sender-specific predictors.
#' @param Xr a n*rr matrix of receiver-specific predictors.
#' @param k an integer value indication dimension of multiplicative effect.
#' @param NS number of MCMC iterations. 
#' @param odens number of iterations between saved samples.
#' @param burn size of burn-in period for MCMC.
#' @param priors list of priors. Default value is NULL.
#' @param startv list of starting values. Default value is NULL.
#' @param pred logical indicating whether to sample from predicted distribution. Default is FALSE.
#' @param rho.calc logical indicating whether to estimate rho. Default is FALSE.
#' @param xInclImpList logical indicating whether a list of imputed covariates are supplied. Default is FALSE.
#' @param Xd_L list of imputed dyadic covariates. 
#' @param Xs_L list of imputed sender-specific covariates.
#' @param Xr_L list of imputed receiver-specific covariates.
#' @export

pgbme <-function(
  y = y,  
  Xd = array(dim=c(nrow(Y.start(y)),nrow(Y.start(y)),0)),
  Xs = matrix(nrow=nrow(Y.start(y)),ncol=0),
  Xr = matrix(nrow=nrow(Y.start(y)),ncol=0), 
  
  #model specification
  k = 1,
  priors = NULL, # the list of priors
  startv = NULL, # the list of starting values
  burn = 0, #the size of burn-in
  pred = FALSE, #sample from predicted distribution; NO by default
  rho.calc = FALSE, # Don't estimate rho by default
  #Details of output
  NS = 100,         #number of scans of mcmc to run
  odens = 1,         #output density
  xInclImpList = FALSE,
  Xd_L = NULL,
  Xs_L = NULL,
  Xr_L = NULL
)
{
       
  # Dimensions of everything
  Y  <- Y.start(y)
  n  <- ncol(Y)
  rd <- dim(Xd)[3]
  rs <- dim(Xs)[2]
  rr <- dim(Xr)[2]   
  
  # Design matrix for unit-specific predictors
  if(!xInclImpList){
    X.u <- cbind(rep(.5, 2*n), adiag(Xs,Xr))
  } else {
    X.u_L = lapply(1:length(Xs_L), function(i){
      cbind(rep(.5,2*n), adiag(Xs_L[[i]], Xr_L[[i]]))  })
    rm(list=c('Xs_L','Xr_L'))
  }
  
  # Construct an upper triangular matrix (useful for later computations)
  tmp<-matrix(1,n,n)
  tmp[lower.tri(tmp,diag=T)]<-NA 
  UT<-tmp

  # set priors
  if (is.null(priors)){
    pi.Sab <- c(1, 0, 1, 4) # inv-wishart
    pi.s2u <- c(1, 1) # inv-gamma: appriori mean approx 1  
    pi.s2v <- c(1, 1) # inv-gamma: appriori mean approx 1
    pi.s2z <- matrix(c(1, 1), max(k, 1), 2, byrow=TRUE) # inv-gamma for z's
    pim.bd <- rep(0, rd) 
    piS.bd <- 10*diag(rd)
    pim.b0sr <- rep(0, 1 + rs + rr) 
    piS.b0sr <- 10*diag(1 + rs + rr) 
  }
  
  # get starting values
  if (is.null(startv)){  
    X <- apply(Xd, 3, c) 
    s <- rep(1:n, n)
    r <- rep(1:n, each=n)           
    cat("Using MLE to calculate starting values", '\n')
    options(warn=-1)    
    mle  <- glmer(c(Y) ~ X + (1|s) + (1|r), 
      family=binomial(link='probit'), nAGQ=0,
      control=glmerControl( optimizer = "nloptwrap" ) )    
    startv$beta.d <- fixef(mle)[-1]
    startv$beta.u <- c(fixef(mle)[1], rep(0, rs+rr))
    startv$s <- ranef(mle)$s[,1]
    startv$r <- ranef(mle)$r[,1]
    startv$z <- matrix(0, n, max(k,1))
  }
  
  # org results
  beta.d<-startv$beta.d ; beta.u<-startv$beta.u 
  s<-startv$s ; r<-startv$r ; z<-startv$z;
  rm(list=c('X','mle','startv'))
  e <- z; f <- z;
  su <- sv <- pi.s2u[2]
  rho <- 0; se <-1;

  ###Matrices for ANOVA  on Xd, s, r
  tmp<-TuTv(ncol(Y))   #unit level effects design matrix for u=yij+yji, v=yij-yji
  Tu<-tmp$Tu
  Tv<-tmp$Tv
  
  #regression design matrix for u=yij+yji, v=yij-yji
  if(!xInclImpList){
    tmp<-XuXv(Xd)   
    Xu<-tmp$Xu  
    Xv<-tmp$Xv
    
    XTu<-cbind(Xu,Tu)
    XTv<-cbind(Xv,Tv)
    
    tXTuXTu<-t(XTu)%*%XTu
    tXTvXTv<-t(XTv)%*%XTv
    rm(list=c('Xu','Xv','tmp'))
  } else {
    tmp <- list()
    for(ii in 1:length(Xd_L)){ tmp[[ii]] <- XuXv(Xd_L[[ii]]) }
    Xu_L <- lapply(tmp, function(x){x$Xu})
    Xv_L <- lapply(tmp, function(x){x$Xv})
    
    XTu_L <- lapply(Xu_L, function(x){cbind(x,Tu)})
    XTv_L <- lapply(Xv_L, function(x){cbind(x,Tv)})
    rm(list=c('Xu_L','Xv_L','tmp'))    

    tXTuXTu_L <- lapply(XTu_L, function(x){ t(x)%*%x })
    tXTvXTv_L <- lapply(XTv_L, function(x){ t(x)%*%x })    
  }
  
  ###redefining hyperparams
  Sab0<-matrix( c(pi.Sab[1],pi.Sab[2],pi.Sab[2],pi.Sab[3]),nrow=2,ncol=2) 
  v0<-pi.Sab[4]   
  
  ###initializing error matrix
  E <- matrix(0,nrow=n, ncol=n)

  theta <- Y.hat <- Y.start(y)
  Y     <- cbind(y, y)
  Z     <- rmnorm(n, varcov = diag(n)) + theta
  
  theta.l <- mat.vect(theta)[,1]
  theta.u <- mat.vect(theta)[,2]
  
  #column names for output file
  cnames <- c("ll",
    paste("bd", 1:rd, sep="-"), "Intercept", 
    paste("bs", 1:rs, sep=""), 
    paste("br", 1:rr, sep=""), 
    "s2a","sab","s2b","se2","rho",
    paste("s2e", 1:k,sep=""),
    paste("s2f", 1:k,sep="")
  )
  
  yhat   <- matrix(NA, ncol = (NS-burn)/odens, nrow = n^2)
  mcmc.samp <- matrix(NA, nrow = (NS-burn)/odens, ncol = length(cnames))
  mcmc.e <- array(NA, dim = c(n, k, (NS-burn)/odens))
  mcmc.f <- array(NA, dim = c(n, k, (NS-burn)/odens))
  mcmc.s <- matrix(NA, nrow = (NS-burn)/odens, ncol = n)
  mcmc.r <- matrix(NA,  nrow = (NS-burn)/odens, ncol = n)
  mcmc.xd <- rep(NA, (NS-burn)/odens) ; samp <- 1
  nst  <- 1
  
  # saving latent space measurements
  uList = list()
  vList = list()
  latSpaceCntr = 1

  main.time <- proc.time()

  ##### The Markov chain Monte Carlo algorithm
  for(ns in 1:NS){  
  
    if(xInclImpList){
      ## sample from posterior of imputed data
      samp <- sample(1:length(X.u_L), 1)
      
      # dyadic and unit level data
      Xd <- Xd_L[[samp]]
      X.u <- X.u_L[[samp]]

      # pre calc'd reg design arrays
      XTu<-XTu_L[[samp]]
      XTv<-XTv_L[[samp]]
      tXTuXTu<-tXTuXTu_L[[samp]]
      tXTvXTv<-tXTvXTv_L[[samp]]
    }

    ###impute any missing values 
    if(any(is.na(Y))) {  
      mu<-theta.betaX.d.srE.ef(beta.d,Xd,s,r,E*0,e,f) #predicted means
      theta.new=mat.vect(theta)
      mu.new=mat.vect(mu)
      t = which(is.na(theta.new[,1]) & is.na(theta.new[,2]))
      theta.new[t,] = rmnorm(length(t), c(0, 0), matrix(se*c(1,rho,rho,1),2,2)) + mu.new[t,]
      t=which(is.na(theta.new[,1]))
      theta.new[t,1] = rnorm(
        length(t), mu.new[t,1] + rho*(mu.new[t,2] - theta.new[t,2]), sqrt(se*(1-rho^2)))
      t=which(is.na(theta.new[,2]))
      theta.new[t,2] = rnorm(
        length(t), mu.new[t,2] + rho*(mu.new[t,1] - theta.new[t,1]), sqrt(se*(1-rho^2)))
      theta.new = vect.mat(theta.new)
      diag(theta.new) = diag(theta)
      theta = theta.new
      
      Y.hat[is.na(Y.start(Y[,1]))] = 1*(theta[is.na(Y.start(Y[,1]))]>0)
      y <- apply(mat.vect(Y.hat), 1, prod)
      Z <- rmnorm(n, varcov = diag(n)) + theta
    }
    ## END IMPUTING MISSING VALUES  
    
    ## AUGMENT LATENT VARIABLES FOR THE PROBIT MODEL
    mu  <- theta.betaX.d.srE.ef(beta.d, Xd, s, r, E*0, e, f)
    Z <- update.Z(y, mu, rho, Z)
    theta <- Z
   
    ### Update regression part
    tmp <- uv.E(theta-e%*%t(f), UT) #the regression part
    u <- tmp$u                  #u=yij+yji,  i<j
    v <- tmp$v                  #v=yij-yji,  i<j

    # Update correlation rho:
    if(rho.calc){
    E <- theta - theta.betaX.d.srE.ef(beta.d, Xd, s, r, E*0, e, f) # update error matrix

      rho.post <- function(z) {
        rho <- fisher.inv(z)
        sum(dmnorm(mat.vect(E), varcov = covmat(rho), log=TRUE)) + dnorm(z, 0, 3/4, log=TRUE)
      }  
      
      cand <- fisher(rho) + rnorm(1, 0, 0.01*log(n))
      a <- min(0,  rho.post(cand) - rho.post(fisher(rho)))
      if (runif(1) < exp(a)) rho <- fisher.inv(cand)
    }
    
    se <- 1
    sv <- 2 - 2*rho
    su <- 4 - sv

    sr.hat <- X.u%*%beta.u    #Sab  
    a <- s-sr.hat[1:n]
    b <- r-sr.hat[n+(1:n)]
    Sab <- rSab.gibbs(a,b,Sab0,v0)
    
    #dyad specific regression coef and unit level effects
    mu<-c(pim.bd, X.u%*%beta.u)    #"prior" mean for (beta.d,s,r)
    beta.d.sr <- rbeta.d.sr.gibbs(u,v,su,sv,piS.bd,mu,Sab,n,XTu,XTv,tXTuXTu,tXTvXTv) 
    beta.d <- beta.d.sr$beta.d ; s<-beta.d.sr$s ; r<-beta.d.sr$r
    #regression coef for unit level effects
    beta.u<-rbeta.sr.gibbs(s,r,X.u,pim.b0sr,piS.b0sr,Sab)           
    
    ### bilinear effects Z
    if(k>0){
      
      #update variance
      s2e<-1/rgamma(k,pi.s2z[,1]+n/2,pi.s2z[,2] + diag(t(e)%*%e)/2)
      s2f<-1/rgamma(k,pi.s2z[,1]+n/2,pi.s2z[,2] + diag(t(f)%*%f)/2)
      
      #Gibbs for zs, using regression
      res<-theta-theta.betaX.d.srE.ef(beta.d,Xd,s,r,0*E,0*e,0*f) 
      s2u.res<-2*se*(1+rho)
      s2v.res<-2*se*(1-rho)
      
      for(i in sample(1:n)){
        u.res<-(res[i,-i]+res[-i,i])
        v.res<-(res[i,-i]-res[-i,i])
        
        alp.j<-cbind( f[-i,],e[-i,] )
        gam.j<-cbind( f[-i,],-e[-i,])
        
        Sef<-chol2inv(
          chol(
            diag(1/c(s2e,s2f),nrow=2*k) + t(alp.j)%*%alp.j/s2u.res + t(gam.j)%*%gam.j/s2v.res ))
        
        muef<-Sef%*%( t(alp.j)%*%u.res/s2u.res + t(gam.j)%*%v.res/s2v.res )
        ef<-t(rmvnorm(muef,Sef))
        e[i,]<-ef[1,1:k] ; f[i,]<-ef[1,k+1:k]
      }
    }
    
    ###
    if(k == 0) {sz <- NULL; s2e <- 0; s2f <- 0}
    
    if(ns > burn & (ns-burn)%%odens==0){  
      n.samp <- (ns-burn)/odens
      E <- theta-theta.betaX.d.srE.ef(beta.d,Xd,s,r,E*0,e, f)
      tmp <- uv.E(theta-e%*%t(f), UT) #the regression part
      u <- tmp$u                  #u=yij+yji,  i<j
      v <- tmp$v                  #v=yij-yji,  i<j
      lpy.th <- sum(dnorm(u,0,sqrt(su),log=T) + dnorm(v,0,sqrt(sv),log=T))
      
      out <- c(lpy.th, beta.d, beta.u, 
        Sab[1,1], Sab[1,2], Sab[2,2], se, rho, s2e[0:k], s2f[0:k])
      
      yhat[, nst] <- c(theta - E)
      mcmc.samp[nst, ] <- t(out)
      mcmc.e[, , nst] <- e
      mcmc.f[, , nst] <- f
      mcmc.s[nst, ] <- s
      mcmc.r[nst, ] <- r
      mcmc.xd[nst] <- samp
      nst         <- nst + 1
      uList[[latSpaceCntr]] <- u
      vList[[latSpaceCntr]] <- v
      latSpaceCntr <- latSpaceCntr + 1

    } # end of output function
    
    if (ns==1){
      cat(paste("MCMC sampling. Estimated time ", 
        format(.POSIXct(
          NS*(proc.time() - main.time)[3],tz="GMT"), "%H:%M:%S"), sep = ""), '\n')
      cat("Progress: ", "\n")
    }
    
    if(ns%%(NS/10)==0){cat(100*ns/NS,"%  ", sep="")}
    
    
  } ## end of MCMC function
  
  cat("\n")
  cat("Time elapsed: ", format(.POSIXct((
    proc.time() - main.time)[3],tz="GMT"), "%H:%M:%S"), sep = " ", "\n")
  #out <- read.table(paste(out.name, "theta", sep ="/"), header = TRUE)
  colnames(mcmc.samp) <- cnames
  mcmc.samp <- list(
    est = mcmc.samp, Xd = Xd, s = mcmc.s, r = mcmc.r, 
    e = mcmc.e, f = mcmc.f, Xd_L=Xd_L, xdId=mcmc.xd,
    uList = uList, vList = vList
    , yhat=yhat
    )
  class(mcmc.samp) <- "gbme"
  return(mcmc.samp)
}