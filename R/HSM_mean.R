# sampler for mean estimation problem
# Gibbs Sampler -----------------------------------------------------------
HS_mix_mean <- function(Y,
                   L = 50,
                   verbose = T,
                   nsim = 3000,
                   burn_in = 3000,
                   thinning = 1,
                   a_dir = .05,
                   BNP = F,
                   a_BNP = 1,
                   d_BNP=1,
                   DP=T,
                   a_a_BNP = 1,
                   b_a_BNP = 1,
                   tau2fixed = 0,
                   gamma.tau=1,
                   gamma.lambda=1,
                   tau_IG=F,
                   cheap = F) {

  n       <- length(Y)

  # Initial values
  a       <- rep(a_dir, L)
  pi      <- MCMCpack::rdirichlet(1, a)
  lambda2 <- rcauchy(L) ^ 2
  zeta    <- sample(L, size = n, replace = T, pi)
  zeta    <- factor(zeta, levels = 1:L)
  lambda2Zeta <- lambda2[zeta]

  nj      <- table(zeta)
  sigma2  <- rgamma(1, 1, 1)

  if (tau2fixed == 0) {
    tau2    <- rgamma(1, 1, 1)
  } else{
    tau2 <- tau2fixed
  }

  random <- F
  if(DP & BNP & a_BNP==0){
    a_BNP  <- rgamma(1,a_a_BNP,b_a_BNP)
    random <- T
  }

  # containers --------------------------------------------------------------
  ZETA    <- array(NA, c(nsim, n))
  BETA    <- array(NA, c(nsim, n))

  PI      <- array(NA, c(nsim, L))
  SIGMA2  <- numeric(nsim)
  TAU2    <- numeric(nsim)
  A_BNP   <- numeric(nsim)

  LAMBDA2 <- array(NA, c(nsim, L))

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


    # Update means
      unosukappa <- 1  + (1 / (lambda2Zeta * tau2 ))
      beta       <- rnorm(n,  (Y ) / (unosukappa), sqrt(sigma2/unosukappa)  )

    # Update sigma2
      e       <- sum((Y - beta)^2)/2
      RSigma  <- .5 * ( sum( (beta ^ 2 )/(lambda2Zeta *tau2) )) + e  #)
      sigma2  <- rivngamma_cpp(a = n, b = RSigma)
    ###################################################################

    if (tau2fixed > 0) {
      tau2 <- tau2fixed
    } else if(tau_IG){
      tau2 <- Tau2_IGprior(beta = beta,lambda2zeta = lambda2Zeta,sigma2 = sigma2,
                           atau =1, btau=1)
    }else{
      tau2    <- SSTau2_gamma2(
        beta = beta,
        lambda2 = lambda2,
        zeta = zeta,
        sigma2 = sigma2,
        oldtau2 = tau2,
        p = n,
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
      zeta <-   update_zeta_cpp(beta = beta,p = n,
                                lambda2 = lambda2,
                                pi = pi,
                                sigma2 = sigma2,
                                L = L,
                                tau2 = tau2,
                                posslab = 1:L)

    zeta    <- factor(zeta, levels = 1:L)
    nj      <- table(zeta)

    if(BNP & DP & random){
      a_BNP <- rgamma(1, a_a_BNP + L - 1, b_a_BNP - sum(log(u[1:(L-1)] )))
    }

    ############################################################################################
    lambda2 <- post.lambda2(lambda2old = lambda2,
                            tau2 = tau2,
                            zeta = zeta,
                            beta = beta,
                            L = L,
                            sigma2 = sigma2,
                            nj = nj)

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
      LAMBDA2[rr,] <- lambda2
    }

    if (verbose) {
      ipbar <- ipbar + 1
      setTxtProgressBar(pbar, ipbar)
    }
  }


  if(cheap){
    out <- list(
      sigma2  = SIGMA2,
      tau2    = TAU2,
      pi      = PI,
      zeta    = ZETA,
      a_BNP   = A_BNP,
      mu_beta    = colMeans(BETA),
      med_beta = apply(BETA,2,median))

    }else{

      out <- list(
        beta    = BETA,
        sigma2  = SIGMA2,
        tau2    = TAU2,
        pi      = PI,
        lambda2 = LAMBDA2,
        zeta    = ZETA,
        a_BNP   = A_BNP,
        mu_beta    = colMeans(BETA),
        Y = Y)


    }



  return(out)
}


