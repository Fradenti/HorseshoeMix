Rcpp::sourceCpp("cpp/HSM.cpp")
source("R/HSM_Auc.R")
source("R/VVS_Gibbs.R")
source('R/Auxiliary_functions/DataGenerator.R')

library(tidyverse)
# A, B, C
# nsim - burnin - thinning - n - p - nb
s1 <- c(5000, 5000, 5, 500, 250, 100)
s2 <- c(5000, 5000, 5, 500, 500, 100)
s3 <- c(5000, 5000, 5, 500, 750, 100)
Scenario <- list(s1,s2,s3)
seedlist <-round( seq(0,5900,length.out=30) )

Y1 <- Data.generator.regression(n = Scenario[[1]][4],
                                p = Scenario[[1]][5],
                                nb1 = Scenario[[1]][6],
                                seedlist = seedlist,
                                NSIM = 30,
                                sdX = 1)
Y2 <- Data.generator.regression(n = Scenario[[2]][4],
                                p = Scenario[[2]][5],
                                nb1 = Scenario[[2]][6],
                                seedlist = seedlist,
                                NSIM = 30,
                                sdX = 1)
Y3 <- Data.generator.regression(n = Scenario[[3]][4],
                                p = Scenario[[3]][5],
                                nb1 = Scenario[[3]][6],
                                seedlist = seedlist,
                                NSIM = 30,
                                sdX = 1)

data_Y <- list(Y1,Y2,Y3)

Dataset_4_Bayesreg_1 <- map(Y1, 
                          ~ as.data.frame(cbind(.x$Y,.x$X)))
Dataset_4_Bayesreg_2 <- map(Y2, 
                            ~ as.data.frame(cbind(.x$Y,.x$X)))
Dataset_4_Bayesreg_3 <- map(Y3, 
                            ~ as.data.frame(cbind(.x$Y,.x$X)))

data_BR <- list(Dataset_4_Bayesreg_1,Dataset_4_Bayesreg_2,Dataset_4_Bayesreg_3)

################################################# parallelizable functions
# HSMix
HSMIX <- function(i,DF_Y){
  HS_mix_regression(Y = DF_Y$Y,
         X = DF_Y$X,
         Intercept = T,
         a_dir = .05,
         cheap = T,
         thinning = 5,
         nsim = nsim,
         burn_in = BI,
         tau2fixed = 0.0001,
         verbose = T)
}
# HS
bayreg_HS <- function(i,DF){
  bayesreg::bayesreg(V1~.,
                     data = DF,model = "normal",
                     prior="hs", n.samples = nsim, burnin = BI,thin = 5)
}
# LASSO
bayreg_LS <- function(i,DF){
  bayesreg::bayesreg(V1~.,
                     data = DF,model = "normal",
                     prior="lasso", n.samples = nsim, burnin = BI,thin = 5)
}
# HS+
bayreg_HSplus <- function(i,DF){
  bayesreg::bayesreg(V1~.,
                     data = DF,model = "normal",
                     prior="hs+", n.samples = nsim, burnin = BI,thin = 5)
}
# RIDGE
bayreg_Ridge <- function(i,DF){
  bayesreg::bayesreg(V1~.,
                     data = DF,model = "normal",
                     prior="ridge", n.samples = nsim, burnin = BI,thin = 5)
}


nsim = 5000
BI = 5000
Res_HS <- Res_HSplus <- Res_Ridge <- Res_Lasso <- list()
jj=1

for(jj in 1:3){
  set.seed(1234*jj)
# Run functions in parallel
  Res_HSmix  <- parallel::mclapply(1:30,FUN = function(x) HSMIX(x,DF_Y = data_Y[[jj]][[x]]), mc.cores = 15, mc.cleanup = T)
  saveRDS(Res_HSmix,paste0("simu_HSmix_",jj,".RDS"))
}

for(jj in 1:3){
    set.seed(1230*jj)
    Res_HS[[jj]]     <- parallel::mclapply(1:30,FUN = function(x) bayreg_HS(x,DF=data_BR[[jj]][[x]]),      mc.cores = 15, mc.cleanup = T)
  }
saveRDS(Res_HS,"Bayereg_HS.RDS")

for(jj in 1:3){
  set.seed(1231*jj)
  Res_Lasso[[jj]]  <- parallel::mclapply(1:30,FUN = function(x) bayreg_LS(x,DF=data_BR[[jj]][[x]]),      mc.cores = 15, mc.cleanup = T)
}
saveRDS(Res_Lasso,"Bayereg_Lasso.RDS")

for(jj in 1:3){
  set.seed(1224*jj)
  Res_HSplus[[jj]] <- parallel::mclapply(1:30,FUN = function(x) bayreg_HSplus(x,DF=data_BR[[jj]][[x]]),      mc.cores = 15, mc.cleanup = T)
}
saveRDS(Res_HSplus,"Bayereg_HSplus.RDS")

for(jj in 1:3){
  set.seed(12132*jj)
  Res_Ridge[[jj]]  <- parallel::mclapply(1:30,FUN = function(x) bayreg_Ridge(x,DF=data_BR[[jj]][[x]]),      mc.cores = 15, mc.cleanup = T)
}
saveRDS(Res_Ridge,"Bayereg_Ridge.RDS")