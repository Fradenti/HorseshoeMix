source("code/R/Auxiliary_functions/DataGenerator.R")
library(tidyverse)
#########
strata_mse <- function(x, y, ind) {
  z <- (x - y) ^ 2
  tapply(z, ind, function(q)
    mean(q))

}
##########
mse <- function(a, b) {
  mean((a - b) ^ 2)
}


D1 <-
  Data.generator.mean(
    n = 500,
    seedlist = 1234 + 1:30,
    nb1 = 50,
    NSIM = 30
  )
D2 <-
  Data.generator.mean(
    n = 500,
    seedlist = 1234 + 1:30,
    nb1 = 100,
    NSIM = 30
  )
D3 <-
  Data.generator.mean(
    n = 1000,
    seedlist = 1234 + 1:30,
    nb1 = 50,
    NSIM = 30
  )
D4 <-
  Data.generator.mean(
    n = 1000,
    seedlist = 1234 + 1:30,
    nb1 = 100,
    NSIM = 30
  )

DATA <- list(D1, D2, D3, D4)

res_extract <- function(num, DATA, nb1, n) {
  # change folder according to your path to saved results
  file1 <- paste0("results_", num, "_hms.RDS")
  file2 <- paste0("results_", num, "_hs.RDS")
  file3 <- paste0("results_", num, "_ssl.RDS")
  sim1 <- readRDS(file1)
  sim2 <- readRDS(file2)
  sim3 <- readRDS(file3)

  D1  <- DATA[[num]]
  ind <- rep(1:3, c(nb1, nb1, n - 2 * nb1))

  overall <- list()
  ###################################
  overall[[1]] <-
    overall_mse_hsm_1 <-
    unlist(map2(D1, sim1, ~ mse(.x$betat, .y[[1]]$mu_beta)))
  overall[[2]] <-
    overall_mse_hsm_2 <-
    unlist(map2(D1, sim1, ~ mse(.x$betat, .y[[2]]$mu_beta)))
  overall[[3]] <-
    overall_mse_hsm_3 <-
    unlist(map2(D1, sim1, ~ mse(.x$betat, .y[[3]]$mu_beta)))
  overall[[4]] <-
    overall_mse_hsm_4 <-
    unlist(map2(D1, sim1, ~ mse(.x$betat, .y[[4]]$mu_beta)))
  overall[[5]] <-
    overall_mse_hsm_5 <-
    unlist(map2(D1, sim1, ~ mse(.x$betat, .y[[5]]$mu_beta)))
  overall[[6]] <-
    overall_mse_hsm_6 <-
    unlist(map2(D1, sim1, ~ mse(.x$betat, .y[[6]]$mu_beta)))
  ###################################
  overall[[7]] <-
    overall_mse_hs_1 <-
    unlist(map2(D1, sim2, ~ mse(.x$betat, .y[[1]]$BetaHat)))
  overall[[8]] <-
    overall_mse_hs_2 <-
    unlist(map2(D1, sim2, ~ mse(.x$betat, .y[[2]]$BetaHat)))
  ###################################
  post_mean_ssl    <- map(sim3,  ~ colMeans(.x$beta[-c(1:5000), ]))
  overall[[9]] <-
    overall_mse_ssl  <-
    unlist(map2(D1, post_mean_ssl, ~ mse(.x$betat, .y)))
  ###################################



  R1 <-
    matrix(unlist(map(overall, ~ c(mean(
      .x
    ), sd(
      .x
    )))), 9, 2, byrow = T)


  overall_strata <- list()

  overall_strata[[1]] <-
    mse_hsm_1 <-
    matrix(unlist(map2(
      D1, sim1, ~ strata_mse(.x$betat, .y[[1]]$mu_beta, ind = ind)
    )), 30, 3, byrow = T)
  overall_strata[[2]] <-
    mse_hsm_2 <-
    matrix(unlist(map2(
      D1, sim1, ~ strata_mse(.x$betat, .y[[2]]$mu_beta, ind = ind)
    )), 30, 3, byrow = T)
  overall_strata[[3]] <-
    mse_hsm_3 <-
    matrix(unlist(map2(
      D1, sim1, ~ strata_mse(.x$betat, .y[[3]]$mu_beta, ind = ind)
    )), 30, 3, byrow = T)
  overall_strata[[4]] <-
    mse_hsm_4 <-
    matrix(unlist(map2(
      D1, sim1, ~ strata_mse(.x$betat, .y[[4]]$mu_beta, ind = ind)
    )), 30, 3, byrow = T)
  overall_strata[[5]] <-
    mse_hsm_5 <-
    matrix(unlist(map2(
      D1, sim1, ~ strata_mse(.x$betat, .y[[5]]$mu_beta, ind = ind)
    )), 30, 3, byrow = T)
  overall_strata[[6]] <-
    mse_hsm_6 <-
    matrix(unlist(map2(
      D1, sim1, ~ strata_mse(.x$betat, .y[[6]]$mu_beta, ind = ind)
    )), 30, 3, byrow = T)
  ###################################
  overall_strata[[7]] <-
    mse_hs_1 <-
    matrix(unlist(map2(
      D1, sim2, ~ strata_mse(.x$betat, .y[[1]]$BetaHat, ind = ind)
    )), 30, 3, byrow = T)
  overall_strata[[8]] <-
    mse_hs_2 <-
    matrix(unlist(map2(
      D1, sim2, ~ strata_mse(.x$betat, .y[[2]]$BetaHat, ind = ind)
    )), 30, 3, byrow = T)
  ###################################
  post_mean_ssl <- map(sim3,  ~ colMeans(.x$beta[-c(1:5000), ]))
  overall_strata[[9]] <-
    mse_ssl       <-
    matrix(unlist(map2(
      D1, post_mean_ssl, ~ strata_mse(.x$betat, .y, ind)
    )), 30, 3, byrow = T)


  R2 <-
    array(unlist(map(
      overall_strata, ~ apply(.x, 2, function(g)
        c(mean(g), sd(g)))
    )), c(2, 3, 9))


  D <- rbind(
    cbind(reshape2::melt(mse_hsm_1), ind = 1),
    cbind(reshape2::melt(mse_hsm_2), ind = 2),
    cbind(reshape2::melt(mse_hsm_3), ind = 3),
    cbind(reshape2::melt(mse_hsm_4), ind = 4),
    cbind(reshape2::melt(mse_hsm_5), ind = 5),
    cbind(reshape2::melt(mse_hsm_6), ind = 6),
    cbind(reshape2::melt(mse_hs_1), ind = 7),
    cbind(reshape2::melt(mse_hs_2), ind = 8),
    cbind(reshape2::melt(mse_ssl),  ind = 9)
  )

  D <- as_tibble(D)

  D <- D %>% mutate(
    competitors = case_when(ind <= 6 ~ 1,
                            ind <= 8 ~ 2,
                            TRUE ~ 3),
    ind = case_when(
      ind == 1 ~ "HSM 1",
      ind == 2 ~ "HSM 2",
      ind == 3 ~ "HSM 3",
      ind == 4 ~ "HSM 4",
      ind == 5 ~ "HSM 5",
      ind == 6 ~ "HSM 6",
      ind == 7 ~ "HS 1",
      ind == 8 ~ "HS 2",
      ind == 9 ~ "SnS"
    ),
    Var2 = case_when(Var2 == 1 ~ "Block 1",
                     Var2 == 2 ~ "Block 2",
                     Var2 == 3 ~ "Block 3")
  )

  Q <- ggplot(D) +
    geom_boxplot(aes(
      x = ind,
      y = value,
      group = ind,
      col = factor(competitors)
    )) +
    facet_wrap( ~ Var2, scales = "free") +
    theme_bw() + coord_flip() + xlab("") + ylab("MSE") +
    theme(text = element_text(size = 16), legend.position = "none")

  LL = list(overall = R1,
            strata = R2,
            plot = Q)
  return(LL)

}
