###Preconditions of the main function "gbme"

##objects with no defaults
#Y : a square matrix (n*n)

##objects with defaults
#Xd : an n*n*rd array representing dyad-specific predictors
#Xs : an n*rs matrix of sender-specific predictors
#Xr : an n*rr matrix of receiver-specific predictors
#xMiss: 
#Xd_L : a list of n*n*rd arrays representing dyad-specific predictors
#Xs_L : a list of n*rs matrices of sender-specific predictors
#Xr_L : a list of  n*rr matrices of receiver-specific predictors

# fam : "gaussian" "binomial", or "poisson"
# if fam="binomial" then an (n*n) matrix N of trials is needed
#
# pi.Sab = (s2a0,sab0,s2b0,v0) : parameters of inverse wishart dist
# pi.s2u = vector of length 2: parameters of an inv-gamma distribution
# pi.s2v = vector of length 2: parameters of an inv-gamma distribution
# pi.s2z = k*2 matrix : parameters of k inv-gamma distributions
#
# pim.bd  : vector of length rd, the prior mean for beta.d
# piS.bd  : a rd*rd pos def matrix, the prior covariance for beta.d
#
# pim.b0sr : vector of length 1+rs+rr, prior mean for (beta0,beta.s,beta.r)
# piS.b0sr  : a (1+rs+rr)*(1+rs+rr) pos def matrix, 
#             the prior variance for (beta0,beta.s,beta.r)


# starting values
# beta.d   #a vector of length rd
# beta.u   #a vector of length 1+rs+rr
# s        #a vector of length n
# r        #a vector of length n 
# z        #an n*k matrix

###the main  function

library(magic)
library(msm)
library(lme4)
library(mnormt)


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
  seed = 1,
  xMiss = FALSE,
  Xd_L = NULL,
  Xs_L = NULL,
  Xr_L = NULL
)
{

  # set seed
  set.seed(seed)
                                                 
  # Dimensions of everything
  Y  <- Y.start(y)
  n  <<- ncol(Y)
  rd <<- dim(Xd)[3]
  rs <<- dim(Xs)[2]
  rr <<- dim(Xr)[2]   
  
  # Design matrix for unit-specific predictors
  if(!xMiss){
    X.u <- cbind(rep(.5, 2*n), adiag(Xs,Xr))
  } else {
    X.u_L = lapply(1:length(Xs_L), function(i){
      cbind(rep(.5,2*n), adiag(Xs_L[[i]], Xr_L[[i]]))  })
    rm(list=c('Xs_L','Xr_L'))
  }
  
  # Construct an upper triangular matrix (useful for later computations)
  tmp<-matrix(1,n,n)
  tmp[lower.tri(tmp,diag=T)]<-NA 
  UT<<-tmp
  
  # get starting values
  if (is.null(startv)){
    Xdyad <- apply(Xd, 3, function(x) c(mat.vect(x)))
    s <- rep(1:n, n-1)
    r <- rep(1:n, each=n-1)
    cat("Using MLE to calculate starting values", '\n')
    options(warn=-1)
    mle  <- glmer(c(mat.vect(Y)) ~ Xdyad + Xs[s, ] + Xr[r, ] + (1|s) + (1|r), 
      family=binomial(link=probit),
      control=glmerControl( optimizer = "nloptwrap" ) )
    startv$beta.d <- fixef(mle)[grepl("Xdyad", names(fixef(mle)))]
    startv$beta.u <- c(
      fixef(mle)[1], fixef(mle)[grepl("Xs", names(fixef(mle)))], 
      fixef(mle)[grepl("Xr", names(fixef(mle)))])
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

  ###Matrices for ANOVA  on Xd, s, r
  tmp<-TuTv(ncol(Y))   #unit level effects design matrix for u=yij+yji, v=yij-yji
  Tu<-tmp$Tu
  Tv<-tmp$Tv
  
  #regression design matrix for u=yij+yji, v=yij-yji
  if(!xMiss){
    tmp<-XuXv(Xd)   
    Xu<-tmp$Xu  
    Xv<-tmp$Xv
    rm(list=c('Xu','Xv','tmp'))    
    
    XTu<<-cbind(Xu,Tu)
    XTv<<-cbind(Xv,Tv)
    
    tXTuXTu<<-t(XTu)%*%XTu
    tXTvXTv<<-t(XTv)%*%XTv
  } else {
    tmp <- list()
    for(ii in 1:length(Xd_L)){
      tmp[[ii]] <- XuXv(Xd_L[[ii]]) }
    Xu_L <- lapply(tmp, function(x){x$Xu})
    Xv_L <- lapply(tmp, function(x){x$Xv})
    rm(list=c('Xu_L','Xv_L','tmp'))        
    
    XTu_L <- lapply(Xu_L, function(x){cbind(x,Tu)})
    XTv_L <- lapply(Xv_L, function(x){cbind(x,Tv)})

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
  
  mcmc.samp <- matrix(NA, nrow = (NS-burn)/odens, ncol = length(cnames))
  mcmc.e <- array(NA, dim = c(n, k, (NS-burn)/odens))
  mcmc.f <- array(NA, dim = c(n, k, (NS-burn)/odens))
  mcmc.s <- matrix(NA, nrow = (NS-burn)/odens, ncol = n)
  mcmc.r <- matrix(NA,  nrow = (NS-burn)/odens, ncol = n)
  mcmc.xd <- rep(NA, (NS-burn)/odens)
  nst  <- 1
  
  main.time <- proc.time()

  ##### The Markov chain Monte Carlo algorithm

  for(ns in 1:NS){  
  
    if(xMiss){
      ## sample from posterior of imputed data
      samp <- sample(1:length(X.u_L), 1)
      
      # dyadic and unit level data
      Xd <- Xd_L[[samp]]
      X.u <- X.u_L[[samp]]

      # pre calc'd reg design arrays
      XTu<<-XTu_L[[samp]]
      XTv<<-XTv_L[[samp]]
      tXTuXTu<<-tXTuXTu_L[[samp]]
      tXTvXTv<<-tXTvXTv_L[[samp]]
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
    tmp <- uv.E(theta-e%*%t(f)) #the regression part
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
    beta.d.sr <- rbeta.d.sr.gibbs(u,v,su,sv,piS.bd,mu,Sab) 
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
      tmp <- uv.E(theta-e%*%t(f)) #the regression part
      u <- tmp$u                  #u=yij+yji,  i<j
      v <- tmp$v                  #v=yij-yji,  i<j
      lpy.th <- sum(dnorm(u,0,sqrt(su),log=T) + dnorm(v,0,sqrt(sv),log=T))
      
      out <- c(lpy.th, beta.d, beta.u, 
        Sab[1,1], Sab[1,2], Sab[2,2], se, rho, s2e[0:k], s2f[0:k])
      
      mcmc.samp[nst, ] <- t(out)
      mcmc.e[, , nst] <- e
      mcmc.f[, , nst] <- f
      mcmc.s[nst, ] <- s
      mcmc.r[nst, ] <- r
      mcmc.xd[nst] <- samp
      nst         <- nst + 1
  
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
    e = mcmc.e, f = mcmc.f, Xd_L=Xd_L, xdId=mcmc.xd)
  class(mcmc.samp) <- "gbme"
  return(mcmc.samp)
} 

# End of MCMC: below are helper functions

# calculate y hat from pgbme output
calc_yhat <- function(m, xMiss=FALSE){
  # Dyadic coefficients:
  b <- m$est[,grep("bd", colnames(m$est))]
  # Empty zero matrix
  E <- matrix(0, nrow = nrow(m$Xd), ncol = nrow(m$Xd))  
  # Calculate predictions and collapse
  if(!xMiss){
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

# SAMPLE Y from posterior predictive distribtuion
pred.y <- function(theta, rho, se, fam="gaussian"){
  e <- rmnorm(n*(n-1)/2, varcov = se*covmat(rho))
  e <- vect.mat(e)
  Y <- theta + e
  diag(Y) <- 0
  if (fam=="binomial") Y = 1*(Y>0)
  Y
}


####
TuTv<-function(n){
  Xu<-Xv<-NULL
  for(i in 1:(n-1)){
    tmp<-tmp<-NULL
    if( i >1 ){ for(j in 1:(i-1)){ tmp<-cbind(tmp,rep(0,n-i)) } }
    tmp<-cbind(tmp,rep(1,n-i)) 
    tmpu<-cbind(tmp,diag(1,n-i)) ; tmpv<-cbind(tmp,-diag(1,n-i))
    Xu<-rbind(Xu,tmpu) ; Xv<-rbind(Xv,tmpv)
  }
  
  list(Tu=cbind(Xu,Xu),Tv=cbind(Xv,-Xv))
}
####

XuXv<-function(X){
  Xu<-Xv<-NULL
  if(dim(X)[3]>0){
    for(r in 1:dim(X)[3]){
      xu<-xv<-NULL
      for(i in 1:(n-1)){
        for(j in (i+1):n){ xu<-c(xu,X[i,j,r]+X[j,i,r])
                           xv<-c(xv,X[i,j,r]-X[j,i,r]) }}
      Xu<-cbind(Xu,xu)
      Xv<-cbind(Xv,xv)  } 
  }
  list(Xu=Xu,Xv=Xv)}

###
uv.E<-function(E){
  u<- c(  t( (  E + t(E))*UT ) )
  u<-u[!is.na(u)]
  v<-c(  t( (  E - t(E))*UT ) )
  v<-v[!is.na(v)]
  list(u=u,v=v)
}
####


####
theta.betaX.d.srE.z<-function(beta.d,X.d,s,r,E,z){
  m<-dim(X.d)[3]
  mu<-matrix(0,nrow=length(s),ncol=length(s))
  if(m>0){for(l in 1:m){ mu<-mu+beta.d[l]*X.d[,,l] }}
  tmp<-mu+re(s,r,E,z)
  diag(tmp)<-0
  tmp
}
####


theta.betaX.d.srE.ef<-function(beta.d,X.d,s,r,E,e,f){
  m<-dim(X.d)[3]
  mu<-matrix(0,nrow=length(s),ncol=length(s))
  if(m>0){for(l in 1:m){ mu<-mu+beta.d[l]*X.d[,,l] }}
  tmp<-mu+reef(s,r,E,e,f)
  diag(tmp)<-0
  tmp}


####
re<-function(a,b,E,z){
  n<-length(a)
  matrix(a,nrow=n,ncol=n,byrow=F)+matrix(b,nrow=n,ncol=n,byrow=T)+E+z%*%t(z) 
}
####

reef<-function(a,b,E,e,f){
  n<-length(a)
  matrix(a,nrow=n,ncol=n,byrow=F)+matrix(b,nrow=n,ncol=n,byrow=T)+E+e%*%t(f) }
####


####
rbeta.d.sr.gibbs <- function(u,v,su,sv,piS.bd,mu,Sab){
  del<-Sab[1,1]*Sab[2,2]-Sab[1,2]^2
  iSab<-rbind(cbind( diag(rep(1,n))*Sab[2,2]/del ,-diag(rep(1,n))*Sab[1,2]/del),
              cbind( -diag(rep(1,n))*Sab[1,2]/del,diag(rep(1,n))*Sab[1,1]/del) )
  rd<-dim(as.matrix(piS.bd))[1]
  
  if(dim(piS.bd)[1]>0){
    cov.beta.sr<-matrix(0,nrow=rd,ncol=2*n)
    iS<-rbind(cbind(solve(piS.bd),cov.beta.sr),
              cbind(t(cov.beta.sr),iSab)) }
  else{iS<-iSab}
  
  Sig <- chol2inv(chol(iS + tXTuXTu/su + tXTvXTv/sv))
  #this may have a closed form expression
  
  M<-Sig%*%(t((u%*%XTu)/su + (v%*%XTv)/sv) + iS%*%mu)
  
  beta.sr<-rmvnorm(M, Sig)
  list(beta.d=beta.sr[(rd>0):rd],s=beta.sr[rd+1:n],r=beta.sr[rd+n+1:n]) }
####


####
rse.beta.d.gibbs<-function(g0,g1,x,XTx,s,r,beta.d){
  n<-length(s)
  1/rgamma(1, g0+choose(n,2)/2,g1+.5*sum( (x-XTx%*%c(beta.d,s,r))^2 ) ) }
####


####
rbeta.sr.gibbs<-function(s,r,X.u,pim.b0sr,piS.b0sr,Sab) {
  del<-Sab[1,1]*Sab[2,2]-Sab[1,2]^2
  iSab<-rbind(cbind( diag(rep(1,n))*Sab[2,2]/del ,-diag(rep(1,n))*Sab[1,2]/del),
              cbind( -diag(rep(1,n))*Sab[1,2]/del,diag(rep(1,n))*Sab[1,1]/del) )
  
  S<-solve( solve(piS.b0sr) + t(X.u)%*%iSab%*%X.u )
  mu<-S%*%(  solve(piS.b0sr)%*%pim.b0sr+ t(X.u)%*%iSab%*%c(s,r))
  rmvnorm( mu,S)
}
####

####
rSab.gibbs<-function(a,b,S0,v0){
  n<-length(a)
  ab<-cbind(a,b)
  Sn<-S0+ (t(ab)%*%ab)
  solve(rwish(solve(Sn),v0+n) )
}

####
rmvnorm<-function(mu,Sig2){
  R<-t(chol(Sig2))
  R%*%(rnorm(length(mu),0,1)) + mu}


####
rwish<-function(S0,nu){ 
  S<-S0*0
  for(i in 1:nu){ z<-rmvnorm(rep(0,dim(as.matrix(S0))[1]), S0)
                  S<-S+z%*%t(z)  }
  S 
}

###  Procrustes transformation: rotation and reflection
proc.rr<-function(Y,X){
  k<-dim(X)[2]
  A<-t(Y)%*%(X%*%t(X))%*%Y
  eA<-eigen(A,symmetric=T)
  Ahalf<-eA$vec[,1:k]%*%diag(sqrt(eA$val[1:k]),nrow=k)%*%t(eA$vec[,1:k])
  t(t(X)%*%Y%*%solve(Ahalf)%*%t(Y)) 
}

####
rbeta.d.s.gibbs <- function(u,su,piS.bd,mu,s2a){
  iSa<-diag(rep(1/s2a,n))
  
  if(dim(piS.bd)[1]>0){
    cov.beta.s<-matrix(0,nrow=rd,ncol=n)
    iS<-rbind(cbind(solve(piS.bd),cov.beta.s),
              cbind(t(cov.beta.s),iSa)) }
  else{iS<-iSa}
  
  Sig<-chol2inv(chol(iS + tXTuXTu/su))
  #this may have a closed form expression
  
  theta <- Sig%*%(t((u%*%XTu)/su) + iS%*%mu)
  
  beta.s <- rmvnorm(theta,Sig)
  list(beta.d=beta.s[(rd>0):rd],s=beta.s[rd+1:n]) }
####

####
rbeta.s.gibbs<-function(s,X.u,pim.b0s,piS.b0s,s2a) {
  iSa<-diag(rep(1/s2a,n))
  S<-solve( solve(piS.b0s) + t(X.u)%*%iSa%*%X.u )
  mu<-S%*%(  solve(piS.b0s)%*%pim.b0s+ t(X.u)%*%iSa%*%s)
  rmvnorm( mu,S)
}



mat.vect <- function(Y){ ## transform NxN matrix into N*(N-1)x2 vector
  Y.l = t(Y)[upper.tri(t(Y))]
  Y.u = Y[upper.tri(Y)]
  Y.new=cbind(Y.l, Y.u)
  Y.new
}

vect.mat <- function(Y){ ## back-transform N*(N-1)x2 vector into NxN matrix
  n = (sqrt((8*nrow(Y)+1))+1)/2
  M = matrix(0, n, n)
  M[t(lower.tri(M))] = Y[,1]
  M = t(M)
  M[upper.tri(M)] = Y[,2] 
  M
}

covmat <- function(rho) matrix(c(1, rho, rho, 1), 2, 2)

Y.start <- function(y){
  Y <- cbind(y, y)
  Y <- vect.mat(Y)
  diag(Y) <- 0
  Y
}

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
  
  Y
}



#Update Z without sampling Y (for partially observed probit):

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
  
  vect.mat(Z)
  
}

fisher <- function(rho) .5*log((1+rho)/(1-rho))
fisher.inv <- function(z) (exp(2*z) - 1)/(exp(2*z) + 1)
