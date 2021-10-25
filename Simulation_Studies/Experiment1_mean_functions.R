# Source Codes and Libraries ----------------------------------------------
Rcpp::sourceCpp("cpp/HSM.cpp")
source("R/Auxiliary_functions/DataGenerator.R")
source("R/HSM_AUX.R")
source("R/HSM_mean.R")

library(parallel)
library(SSLASSO)
library(horseshoe)
library(BoomSpikeSlab)
library(tidyverse)
library(EnvStats)
library(broom)
library(multcomp)
library(formattable)
library(knitr)
library(ggrepel)
library(GGally)
library(Rtsne)
library(patchwork)
library(latex2exp)


mse <- function(a,b){
  mean((a-b)^2)
}

partitioning_mse <- function(a,b){
  true_big_beta    <- a[1:nb1]
  true_middle_beta <- a[1:nb1+nb1]
  true_small_beta  <- a[-c(1:(2*nb1))]
  est_big_beta     <- b[1:nb1]
  est_middle_beta  <- b[1:nb1+nb1]
  est_small_beta   <- b[-c(1:(2*nb1))]
  
  c(mse(true_big_beta,est_big_beta),
    mse(true_middle_beta,est_middle_beta),
    mse(true_small_beta,est_small_beta))
}


# Simulation Study on Mean ------------------------------------------------

estimate_mean_HSMix <- function(i,D){
  
set.seed(12345*i)
# fixed tau
hsmix_1=HS_mix_mean(Y = D[[i]]$Y,
                     L = 50,
                     a_dir = 0.05,
                     nsim = 5000,
                     burn_in = 5000,
                     tau2fixed = 0.001  ,
                     cheap = T)
# tau IG
hsmix_2=HS_mix_mean(Y = D[[i]]$Y,
                     L = 50,
                     a_dir = 0.05,
                     nsim = 5000,
                     burn_in = 5000,
                     tau2fixed = 0,
                     tau_IG = T  ,
                     cheap = T)
# tau hC
hsmix_3=HS_mix_mean(Y = D[[i]]$Y,
                    L = 50,
                    a_dir = 0.05,
                    nsim = 5000,
                    burn_in = 5000,
                    tau2fixed = 0,
                    tau_IG = F,
                    cheap = T)
# BNP, fixed tau
hsmix_4=HS_mix_mean(Y = D[[i]]$Y,
                    L = 50,
                    BNP = T,a_BNP = 0,
                    a_dir = 0.05,
                    nsim = 5000,
                    burn_in = 5000,
                    tau2fixed = 0.001  ,
                    cheap = T)
# BNP, tau IG
hsmix_5=HS_mix_mean(Y = D[[i]]$Y,
                    L = 50,
                    BNP = T,a_BNP = 0,
                    a_dir = 0.05,
                    nsim = 5000,
                    burn_in = 5000,
                    tau2fixed = 0,
                    tau_IG = T  ,
                    cheap = T)
# BNP, tau hC
hsmix_6=HS_mix_mean(Y = D[[i]]$Y,
                    L = 50,
                    BNP = T,a_BNP = 0,
                    a_dir = 0.05,
                    nsim = 5000,
                    burn_in = 5000,
                    tau2fixed = 0,
                    tau_IG = F  ,
                    cheap = T)
list(hsmix_1,hsmix_2,hsmix_3,hsmix_4,hsmix_5,hsmix_6)
}

estimate_mean_HSpackage <- function(i,D){
  
hs_tau_hC <- horseshoe::HS.normal.means(D[[i]]$Y,                                         
                                        method.sigma = "Jeffreys",
                                        method.tau = "halfCauchy",
                                        burn = 5000,
                                        nmc = 5000)
hs_tau_fix <- horseshoe::HS.normal.means(D[[i]]$Y,
                                         method.sigma = "Jeffreys",
                                         method.tau = "fixed",
                                         burn = 5000,
                                         nmc = 5000,
                                         tau = sqrt(0.0001))

list(hs_tau_hC,hs_tau_fix)
}

estimate_mean_SSL <- function(i,D){
  
ssl <- BoomSpikeSlab::lm.spike(D[[i]]$Y~ - 1 + diag(length(D[[i]]$Y)),
                               niter = 10000)

ssl
}
