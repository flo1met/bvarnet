# Integration test: sim_bvarnet() output → to_stan_data() pipeline
source("R/sim_bvarnet.R")
source("R/to_stan_data.R")

cat("=== Integration 1: Bernoulli → to_stan_data ===\n")
sim <- sim_bvarnet(N = 15, T_obs = 50, p = 3, q = 2, K = 1,
                   family = "bernoulli", seed = 42)
sd <- to_stan_data(data = sim$data, family = "bernoulli",
                   id_col = "id", time_col = "t",
                   y_cols = paste0("y_", 1:3), x_cols = paste0("x_", 1:2),
                   K = 1)
cat("  p:", sd$p, " J:", sd$J, " K:", sd$K, " n_obs:", sd$n_obs, "\n")
cat("  n_fe:", sd$n_fe, " n_re:", sd$n_re, "\n")
cat("  X dim:", dim(sd$X), " B dim:", dim(sd$B), "\n")
cat("  Y range:", range(sd$Y), "\n\n")

cat("=== Integration 2: Gaussian K=2 → to_stan_data ===\n")
sim2 <- sim_bvarnet(N = 10, T_obs = 40, p = 2, q = 1, K = 2,
                    family = "gaussian", seed = 99)
sd2 <- to_stan_data(data = sim2$data, family = "gaussian",
                    id_col = "id", time_col = "t",
                    y_cols = paste0("y_", 1:2), x_cols = "x_1",
                    K = 2)
cat("  p:", sd2$p, " J:", sd2$J, " K:", sd2$K, " n_obs:", sd2$n_obs, "\n")
cat("  n_fe:", sd2$n_fe, " B cols:", ncol(sd2$B), " (expect p*K=4)\n")
cat("  Y dim:", dim(sd2$Y), "\n\n")

cat("=== Integration 3: Ordinal → to_stan_data ===\n")
sim3 <- sim_bvarnet(N = 20, T_obs = 30, p = 2, q = 1, K = 1,
                    family = "ordinal", C = 4, seed = 77)
sd3 <- to_stan_data(data = sim3$data, family = "ordinal",
                    id_col = "id", time_col = "t",
                    y_cols = paste0("y_", 1:2), x_cols = "x_1",
                    K = 1)
cat("  p:", sd3$p, " J:", sd3$J, " C:", sd3$C, "\n")
cat("  n_fe:", sd3$n_fe, " (expect 1: no intercept for ordinal)\n")
cat("  Y range:", range(sd3$Y), "\n\n")

cat("=== Integration 4: re_temporal → to_stan_data ===\n")
sim4 <- sim_bvarnet(N = 10, T_obs = 40, p = 2, q = 1, K = 1,
                    family = "bernoulli", re_temporal = TRUE,
                    sd_phi = 0.2, seed = 5)
sd4 <- to_stan_data(data = sim4$data, family = "bernoulli",
                    id_col = "id", time_col = "t",
                    y_cols = paste0("y_", 1:2), x_cols = "x_1",
                    K = 1, re_temporal = TRUE)
cat("  n_re:", sd4$n_re, " (expect p*K=2 for re_temporal)\n")
cat("  Z dim:", dim(sd4$Z), "\n\n")

cat("=== All integration tests passed! ===\n")
