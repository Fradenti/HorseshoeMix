# Gibbs Sampler -----------------------------------------------------------
#' Horseshoe mix Gibbs sampler
#'
#' @param Y dataset
#' @param X covariates
#' @param L number of mixture components
#' @param Intercept logical, should the intercept be included?
#' @param verbose logical, should I print chain progression?
#' @param cheap logical, should I only save some summaries of the chains
#' @param nsim number of iterations
#' @param burn_in burn in period
#' @param thinning potential thinning period
#' @param a_dir 
#' @param BNP logical, use a nonparametric mixture
#' @param a_BNP hyperparameter for the nonparametric process (DP vs PY)
#' @param d_BNP hyperparameter for the nonparametric process (DP vs PY)
#' @param DP logical, should I use a DP or PY?
#' @param a_a_BNP hyperprior for a_DP \sim Gamma(a_a_BNP,b_a_BNP)
#' @param b_a_BNP hyperprior for a_DP \sim Gamma(a_a_BNP,b_a_BNP)
#' @param tau2fixed logical, should I keep the global shrinkage parameter fixed?
#' @param gamma.tau scale parameters for half Cauchy
#' @param gamma.lambda scale parameters for half Cauchy
#' @param MEAN logical, run faster algorithm for mean estimation
#' @param MEAN.S1 logical, fix sigma2 to 1
#'
#' @return
#' @export
#'
#' @examples
HS_mix_regression <- function(Y,
                   X = NULL,
                   L = 50,
                   Intercept = T,
                   verbose = T,
                   cheap = F,
                   nsim = 3000,
                   burn_in = 3000,
                   thinning = 1,
                   a_dir = 1,
                   BNP = F,
                   a_BNP = 1,
                   d_BNP=1, 
                   DP=T,
                   a_a_BNP = 1,
                   b_a_BNP = 1,
                   tau2fixed = 0,
                   gamma.tau=1, 
                   gamma.lambda=1,
                   MEAN=F, MEAN.S1=F) {
  
  n       <- length(Y)
  if(is.null(X)){ X <- diag(n) }
  p       <- ncol(X)
  
  # precomputables
  if(MEAN==F){
    XTX     <- t(X) %*% (X)
    
    MEAN.S1 <- F
    
    if (p == 1)  XTX     <- c(XTX)
    XTY     <- t(X) %*% (Y)
  }
  
  # Initial values
  a       <- rep(a_dir, L)
  pi      <- MCMCpack::rdirichlet(1, a)
  lambda2 <- rcauchy(L) ^ 2
  zeta    <- sample(L, size = p, replace = T, pi)
  zeta    <- factor(zeta, levels = 1:L)
  lambda2Zeta <- lambda2[zeta]
  
  nj      <- table(zeta)
  sigma2  <- rgamma(1, 1, 1)
  
  if (tau2fixed == 0) {
    tau2    <- rgamma(1, 1, 1)
  } else{
    tau2 <- tau2fixed
  }
  
  if(MEAN){
    c1        <- rep(1,p)
    Intercept <- F
    X         <- diag(n)
    p         <- n
    if(MEAN.S1) {sigma2 <- 1}
  }
  
  if (Intercept) {
    beta0 <- rnorm(1, 0, 10)
  } else{
    beta0 <- 0
  }
  random <- F
  if(DP & BNP & a_BNP==0){
    a_BNP  <- rgamma(1,a_a_BNP,b_a_BNP)
    random <- T
  }
  
  # containers --------------------------------------------------------------
  BETA    <- array(NA, c(nsim, p))
  ZETA    <- array(NA, c(nsim, p))
  
  PI      <- array(NA, c(nsim, L))
  SIGMA2  <- numeric(nsim)
  BETA0   <- numeric(nsim)
  TAU2    <- numeric(nsim)
  A_BNP   <- numeric(nsim)
  
  if (cheap == F) {
    LAMBDA2 <- array(NA, c(nsim, L))
  }
  RUE     <- p / n <= 2
  
  if (verbose) {
    cat("MCMC progress:\n")
    flush.console()
    pbar <- txtProgressBar(min = 1,
                           max = nsim * thinning + burn_in,
                           style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }
  # loop --------------------------------------------------------------------
  for (sim in 1:(nsim * thinning + burn_in)) {
    
    
    if(MEAN){
      unosukappa <- 1  + (1 / (lambda2Zeta * tau2 ))
      beta <- rnorm(p,  (Y - beta0) / (unosukappa), sqrt(sigma2/unosukappa)  )
    }else if (RUE) {
      A       <- XTX / sigma2 
      diag(A) <- diag(A) + (1 / (lambda2Zeta * tau2 * sigma2))
      beta <-
        Sample.MVN.precision.R(Omega = A, 
                                 Xt.Sigma.Y = t(X) %*% (Y - beta0) / sigma2,
                                 p =  p)
      
    } else{
      # BHATTAC. 2016
      beta  <-
        Sample.MVN.precision.A(
          Phi  = X / sqrt(sigma2),
          p = p,
          n = n,
          diagon = (lambda2Zeta * tau2 *
                      sigma2),
          alpha = (Y - beta0) / sqrt(sigma2)
        )
    }
    #########################################################################################
    if (Intercept) {
      beta0 <- rnorm(1, 
                     sum(Y - X %*% beta) / n, 
                     sqrt(sigma2 / n ))
    }
    ###################################################################
    
    if(MEAN.S1){
      sigma2 <- 1
    }else{
      e      <- Y - X %*% beta - beta0
      RSigma <- .5 * ( sum( beta ^ 2 / (lambda2Zeta * (tau2)) ) + sum(e ^2) ) #)
      sigma2  <- 1 / rgamma(1, (n + p) / 2, rate = RSigma)
    }
    ###################################################################
    if (tau2fixed > 0) {
      tau2 <- tau2fixed
    } else{
      tau2    <- SSTau2_gamma2(
        beta = beta,
        lambda2 = lambda2,
        zeta = zeta,
        sigma2 = sigma2,
        oldtau2 = tau2,
        p = p,
        gamma2 = gamma.tau
      )
    }
    ############################################################################################
    if (BNP) {
      if(DP){
        u <- UPD_Sticks_DP_cpp(zeta,L,a_BNP)
        u[L] <- 1  #truncated version
      }else{
        u <- UPD_Sticks_PY_cpp(zeta,L,dPY = d_BNP,alphaPY = a_BNP)
        u[L] <- 1  #truncated version
      }
      pi <- StickBreaker_cpp(u)
    } else{
      pi <- c(MCMCpack::rdirichlet(1, alpha = a + nj))
    }
    
    ############################################################################################
    
    zeta    <- update_zeta_cpp(beta = beta,
                               lambda2 = lambda2,
                               pi = pi,
                               sigma2 = sigma2,
                               L = L,
                               p = p,
                               tau2 = tau2,
                               posslab = 1:L)
    
    zeta    <- factor(zeta, levels = 1:L)
    nj      <- table(zeta)
    
    if(BNP & DP & random){
      a_BNP <- rgamma(1, a_a_BNP + L - 1, b_a_BNP - sum(log(u[1:(L-1)] )))
    }
    
    ############################################################################################
    lambda2 <- Slice_lambda_scaled_cpp(beta=beta,
                                       prevlambda2 = lambda2,
                                       zeta=zeta,nl = nj,
                                       sigma2 = sigma2,
                                       tau2 = tau2,
                                       gamma2_lambda = gamma.lambda,
                                       L = L)
    
    lambda2Zeta <- lambda2[zeta]
    ##################################################
    if (sim > burn_in && ((sim - burn_in) %% thinning == 0)) {
      rr                <- floor((sim - burn_in) / thinning)
      
      BETA[rr,]    <- beta
      SIGMA2[rr]   <- sigma2
      TAU2[rr]     <- tau2
      PI[rr,]      <- pi
      A_BNP[rr]    <- a_BNP
      ZETA[rr,]    <- zeta
      
      if (Intercept) {
        BETA0[rr]  <- beta0
      }
      
      if (cheap == F) {
        LAMBDA2[rr,] <- lambda2
      }
    }
    
    if (verbose) {
      ipbar <- ipbar + 1
      setTxtProgressBar(pbar, ipbar)
    }
  }
  
  if (cheap == F) {
    out <- list(
      beta    = BETA,
      sigma2  = SIGMA2,
      tau2    = TAU2,
      beta0   = BETA0,
      pi      = PI,
      lambda2 = LAMBDA2,
      zeta    = ZETA,
      a_BNP   = A_BNP,
      mu_beta    = colMeans(BETA),
      mu_beta0   = mean(BETA0),
      Y = Y,
      Xstd = X
    )
  } else{
    out <- list(
      beta    = colMeans(BETA),
      beta0   = mean(BETA0),
      sigma2  = SIGMA2,
      tau2    = TAU2,
      pi      = PI,
      zeta    = ZETA,
      a_BNP   = A_BNP,
      Y = Y,
      Xstd = X
    )
  }
  return(out)
}
