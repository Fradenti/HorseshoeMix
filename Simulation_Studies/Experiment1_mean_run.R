# Run the experiment ----------------------------------------------------------

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

RES1_hms = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_HSMix(i, D = D1),
  mc.cores = 15
)
RES1_hs  = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_HSpackage(i, D = D1),
  mc.cores = 15
)
RES1_ssl = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_SSL(i, D = D1),
  mc.cores = 15
)
saveRDS(RES1_hms, "results_1_hms.RDS")
saveRDS(RES1_hs , "results_1_hs.RDS")
saveRDS(RES1_ssl, "results_1_ssl.RDS")
RES2_hms = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_HSMix(i, D = D2),
  mc.cores = 15
)
RES2_hs  = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_HSpackage(i, D = D2),
  mc.cores = 15
)
RES2_ssl = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_SSL(i, D = D2),
  mc.cores = 15
)
saveRDS(RES2_hms, "results_2_hms.RDS")
saveRDS(RES2_hs , "results_2_hs.RDS")
saveRDS(RES2_ssl, "results_2_ssl.RDS")
RES3_hms = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_HSMix(i, D = D3),
  mc.cores = 15
)
RES3_hs  = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_HSpackage(i, D = D3),
  mc.cores = 15
)
RES3_ssl = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_SSL(i, D = D3),
  mc.cores = 15
)
saveRDS(RES3_hms, "results_3_hms.RDS")
saveRDS(RES3_hs , "results_3_hs.RDS")
saveRDS(RES3_ssl, "results_3_ssl.RDS")
RES4_hms = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_HSMix(i, D = D4),
  mc.cores = 15
)
RES4_hs  = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_HSpackage(i, D = D4),
  mc.cores = 15
)
RES4_ssl = parallel::mclapply(
  X = 1:30,
  FUN = function(i)
    estimate_mean_SSL(i, D = D4),
  mc.cores = 15
)
saveRDS(RES4_hms, "results_4_hms.RDS")
saveRDS(RES4_hs , "results_4_hs.RDS")
saveRDS(RES4_ssl, "results_4_ssl.RDS")
