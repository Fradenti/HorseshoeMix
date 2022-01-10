
# VVS functions -----------------------------------------------------------
#######################################################################
# Sample mvn first method ---------------------------------------------------
Sample.MVN.precision.R <- function(Omega,Xt.Sigma.Y,p){
  Lc   = chol(Omega)
  v    = forwardsolve(t(Lc), (Xt.Sigma.Y))
  r1   = backsolve(Lc, v)
  w    = backsolve(Lc, stats::rnorm(p, 0, 1))
  beta = r1 + w
  return(beta)
}

# Sample mvn second method ---------------------------------------------------
Sample.MVN.precision.A <- function(Phi,diagon,alpha,p,n){
  u      = as.matrix(stats::rnorm(p, 0, 1)) * sqrt(diagon)
  delta  = as.matrix(stats::rnorm(n, 0, 1))
  v      = Phi %*% u + delta
  DtPhi  = diag(diagon) %*% t(Phi)
  W      = Phi %*% DtPhi + diag(1, n)
  w      = solve( W, (alpha - v))
  theta  = u + DtPhi %*% w
  return(theta)
}
###########################
Sample.MVN.precision.A2 <- function(Phi,diagon,alpha,p,n){
  u      = as.matrix(stats::rnorm(p, 0, 1)) * sqrt(diagon)
  delta  = as.matrix(stats::rnorm(n, 0, 1))
  v      = Phi %*% u + delta
  Dpt = Phi
  for (i in 1:length(diagon)) {
    Dpt[, i] = Dpt[, i] * diagon[i]
  }
  DtPhi = t(Dpt)
  W      = Phi %*% DtPhi + diag(1, n)
  L = chol(W)
  vv = forwardsolve(t(L), (alpha - v))
  w = backsolve(L, vv)
  theta  = u + DtPhi %*% w
  return(theta)
}

#######################################################################
# Slice Sampler TAU -------------------------------------------------------
SSTau2 <- function(beta, lambda2, zeta, sigma2, oldtau2, p){
  et     = 1/oldtau2
  utau   = stats::runif(1, 0, 1/(1 + et))
  ###
  ubt    = (1 - utau)/utau
  tempt  = .5 * sum( beta^2/(lambda2[zeta]*sigma2) )
  Fubt   = stats::pgamma(ubt, (p + 1)/2, rate = tempt)
  Fubt   = max(Fubt, 1e-08)
  ut     = stats::runif(1, 0, Fubt)
  et     = stats::qgamma(ut,  (p + 1)/2, rate = tempt)
  tau2   = 1/(et)
  return(tau2)
}

#############################################
SSTau2_gamma2 <- function(beta, lambda2, zeta, sigma2, oldtau2, p, gamma2){
  et     = gamma2/oldtau2
  utau   = stats::runif(1, 0, 1/(1 + et))
  ###
  ubt    = (1 - utau)/utau
  tempt  = .5/gamma2 * sum( beta^2/(lambda2[zeta]*sigma2) )
  Fubt   = stats::pgamma(ubt, (p + 1)/2, rate = tempt)
  Fubt   = max(Fubt, 1e-08)
  ut     = stats::runif(1, 0, Fubt)
  et     = stats::qgamma(ut,  (p + 1)/2, rate = tempt)
  tau2   = gamma2/(et)
  return(tau2)
}

#######################################################################
# Slice Sampler TAU -------------------------------------------------------
SSTau2.mvt <- function(beta, lambda2, zeta, sigma2, oldtau2, PJ){
  et     = 1/oldtau2
  utau   = stats::runif(1, 0, 1/(1 + et))
  ###
  ubt    = (1 - utau)/utau
  tempt  = .5 * sum( c(beta^2/sigma2) /(lambda2[zeta]) )
  Fubt   = stats::pgamma(ubt, (PJ + 1)/2, rate = tempt)
  Fubt   = max(Fubt, 1e-08)
  ut     = stats::runif(1, 0, Fubt)
  et     = stats::qgamma(ut,  (PJ + 1)/2, rate = tempt)
  tau2   = 1/(et)
  return(tau2)
}

#######################################################################
# Slice Sampler Lambda -------------------------------------------------------
SSLambda2 <- function(beta, oldlambda2, zeta, sigma2, tau2, p, nj, L){

  t      = 1/(oldlambda2)
  upsi   = stats::runif(L, 0, 1/(1 + t))

  lam1   = tapply(beta, zeta, function(x) sum(x^2)/(2*sigma2*tau2))
  lam1[is.na(lam1)] = 0
  ######################################################################
  UPSI <- (1 - upsi)/upsi
  ######################################################################
  Fub   = stats::pgamma(UPSI, (nj + 1)/2, rate = lam1)
  Fub[Fub < (1e-04)] = 1e-04
  ut = stats::runif(L, 0, Fub)
  et = stats::qgamma(ut, (nj + 1)/2, rate = lam1)
  lambda2 = 1/(et)
  K <- sum(lambda2==Inf)
  lambda2[lambda2==Inf] <- (rcauchy(K))^2
  return(lambda2)
}
######################################################################
SSLambda2_v2 <- function(beta, oldlambda2, zeta, sigma2, tau2, p, nj, L){
  lambda2 = numeric(L)
  for(l in 1:L){
    if(nj[l]==0){
      lambda2[l] <- (rcauchy(1))^2
    }else{
      t    <-  1/(oldlambda2[l])
      upsi <-  stats::runif(1, 0, 1/(1 + t))
      lam1 <-  sum(beta[zeta==l]^2)/(2*sigma2*tau2)
      UPSI <- (1 - upsi)/upsi
      Fub  <-  stats::pgamma(UPSI, (nj[l] + 1)/2, rate = lam1)
      Fub[Fub < (1e-04)] <- 1e-04
      ut <- stats::runif(1, 0, Fub)
      et <- stats::qgamma(ut, (nj[l] + 1)/2, rate = lam1)
      lambda2[l] = 1/(et)
    }
  }
  return(lambda2)
}

#######################################################################
#######################################################################
#######################################################################
#######################################################################

# Evaluate the performance vs ground truth --------------------------------
rsse <- function(x,y){
  sqrt(sum((x-y)^2))
}
sse <- function(x,y){
  (sum((x-y)^2))
}
Asse <- function(x,y){
  (sum(abs(x-y)))
}
#######################################################################
#######################################################################
#######################################################################
#######################################################################
##########################################################################

# DP stick breaking posterior ---------------------------------------------
Stick_breaker <- function(u,L){
  omega      <- u# set the length
  u[L]       <- 1
  u.         <- cumprod(1-u)
  a          <- u[2:L]*u.[1:(L-1)] # omega in the stick-breaking representation
  omega[2:L] <- a
  return(omega)
}
Update.stickbreaking <- function(nj,a_BNP,L){
  Nj    <- sum(nj)-cumsum(nj)
  u     <- apply(cbind(nj,Nj), 1, function(x) rbeta(1, 1+x[1], a_BNP+x[2] )   )
  omega <- Stick_breaker(u,L = L)
  return(cbind(u- 1e-5,omega))
}


#######################################################################

# Postprocessing --------------------------------------------

  # Cutting the estimated dendrogram

  est_clust_hier <- function(PSM,K=2){
    dd <- as.dist(1-PSM)
    hc <- hclust(dd)
    cutree(hc,K)
  }

# compute kappa_j

kapper <- function(m1){
  nsim  <- nrow(m1$beta)
  L     <- trace_back(m1)
  kappa <- L
  for (i in 1:nsim){
    kappa[i,] <- 1/(1+m1$tau2[i] * L[i,])
  }
  return(kappa)
}

# Trace Lambdas disentangling the label switching

trace_back <- function(m1){
  nsim= nrow(m1$lambda2)
  N   = ncol(m1$beta)
  r <- matrix(NA_real_, nsim,N)
  for(i in 1:nsim){
    r[i,] <- m1$lambda2[i,m1$zeta[i,]]
  }
  return(r)
}

# clustering, PSM

clustering <- function(m1,type=1, greedy=F){

  PSM <- mcclust::comp.psm(m1$zeta)
  if(type==1){
    if(greedy){
      CL <- mcclust.ext::minVI(PSM,method = "greedy")$cl
    }else{
      CL <- mcclust.ext::minVI(PSM)$cl
    }
  }else if(type==2){
    CL <- GreedyEPL::MinimiseEPL(m1$zeta)$decision
  }else if(type==3){
    CL <- mcclust.ext::minbinder.ext(PSM)$cl
  }

  return(list(PSM=PSM, CL=CL))
}

# Reset
reset <- function(x){
  as.numeric(as.factor(x))
}

# Evaluate the performance vs ground truth --------------------------------
rsse <- function(x,y){
  sqrt(sum((x-y)^2))
}
sse <- function(x,y){
  (sum((x-y)^2))
}
Asse <- function(x,y){
  (sum(abs(x-y)))
}




# Auxiliary Functions ----------------------------------------------------------

post.lambda2 <- function(lambda2old,tau2,zeta,beta,L,sigma2,nj){

  newl2 <- numeric(L)
  t     <- 1/lambda2old
  u     <- runif(L,0,1/(1+t))
  upp   <- 1/u-1
  for(l in 1:L){

    if(nj[l]==0){
      newl2[l] <- rcauchy(1)^2
    }else{
      B2_l  <- sum(beta[zeta==l]^2)/(2*sigma2*tau2)

      upp_p <- pgamma(upp[l],
                      shape = (nj[l]+1)/2,
                      rate =  B2_l)


      upp_p <- ifelse(upp_p<1e-4,1e-4,upp_p)
      pp    <- runif(1,0,upp_p)
      t__   <- qgamma(pp,shape = (nj[l]+1)/2,rate = B2_l,log.p = 0)
      newl2[l] <- 1/t__
    }
  }
  return(newl2)
}



Tau2_IGprior <- function(beta,
                         sigma2,
                         lambda2zeta,
                         atau=atau,
                         btau=btau){

  b.post <-  sum(beta^2/(lambda2zeta*sigma2))/2 + btau
  a.post <-  length(beta)/2 + atau
  tau2   <-  1/rgamma(1, a.post, 1/b.post)

  return(tau2)
}
