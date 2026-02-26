source("R/sim_bvarnet.R")

cat("=== Test 1: Bernoulli (p=4, q=2, K=1, N=25, T=100) ===\n")
sim_bin <- sim_bvarnet(N = 25, T_obs = 100, p = 4, q = 2, K = 1,
                       family = "bernoulli", seed = 42)
cat("  data dim:", dim(sim_bin$data), "\n")
cat("  columns:", paste(names(sim_bin$data), collapse=", "), "\n")
cat("  Y range:", range(sim_bin$data$y_1), "\n")
cat("  alpha:", round(sim_bin$truth$alpha, 3), "\n")
cat("  Phi diag:", round(diag(sim_bin$truth$Phi), 3), "\n")
cat("  VAR stable:", check_var_stability(sim_bin$truth$Phi, 4, 1), "\n\n")

cat("=== Test 2: Gaussian (p=3, q=1, K=2, N=20, T=80) ===\n")
sim_gauss <- sim_bvarnet(N = 20, T_obs = 80, p = 3, q = 1, K = 2,
                         family = "gaussian", seed = 123)
cat("  data dim:", dim(sim_gauss$data), "\n")
cat("  Y range:", round(range(sim_gauss$data$y_1), 2), "\n")
cat("  sigma:", round(sim_gauss$truth$sigma, 3), "\n")
cat("  Phi dim:", dim(sim_gauss$truth$Phi), "\n")
cat("  VAR stable:", check_var_stability(sim_gauss$truth$Phi, 3, 2), "\n\n")

cat("=== Test 3: Ordinal (p=3, q=0, K=1, N=30, T=50, C=5) ===\n")
sim_ord <- sim_bvarnet(N = 30, T_obs = 50, p = 3, q = 0, K = 1,
                       family = "ordinal", C = 5, seed = 999)
cat("  data dim:", dim(sim_ord$data), "\n")
cat("  Y range:", range(sim_ord$data$y_1), "\n")
cat("  Y table node 1:\n")
print(table(sim_ord$data$y_1))
cat("  kappa node 1:", round(sim_ord$truth$kappa[[1]], 3), "\n\n")

cat("=== Test 4: K=3, with random temporal effects ===\n")
sim_re <- sim_bvarnet(N = 15, T_obs = 60, p = 2, q = 1, K = 3,
                      family = "bernoulli", re_temporal = TRUE,
                      sd_phi = 0.15, seed = 7)
cat("  data dim:", dim(sim_re$data), "\n")
cat("  Phi dim:", dim(sim_re$truth$Phi), "\n")
cat("  Phi_i dim:", dim(sim_re$truth$Phi_i), "\n")
cat("  sd_u dim:", dim(sim_re$truth$sd_u), "\n\n")

cat("=== Test 5: Random covariate slopes ===\n")
sim_rs <- sim_bvarnet(N = 20, T_obs = 80, p = 3, q = 2, K = 1,
                      family = "gaussian", sd_gamma = 0.3, seed = 55)
cat("  sd_u dim:", dim(sim_rs$truth$sd_u), "\n")
cat("  gamma_i dim:", dim(sim_rs$truth$gamma_i), "\n\n")

cat("=== Test 6: Ordinal C=2 (binary-like) ===\n")
sim_bin2 <- sim_bvarnet(N = 10, T_obs = 40, p = 2, K = 1,
                        family = "ordinal", C = 2, seed = 11)
cat("  Y range:", range(sim_bin2$data$y_1), "\n")
print(table(sim_bin2$data$y_1))
cat("  kappa node 1:", round(sim_bin2$truth$kappa[[1]], 3), "\n\n")

cat("=== Test 7: q=0, no covariates ===\n")
sim_nocov <- sim_bvarnet(N = 10, T_obs = 30, p = 2, K = 1,
                         family = "bernoulli", q = 0, seed = 1)
cat("  columns:", paste(names(sim_nocov$data), collapse=", "), "\n\n")

cat("=== Test 8: p=1 univariate ===\n")
sim_uni <- sim_bvarnet(N = 10, T_obs = 50, p = 1, K = 1,
                       family = "gaussian", q = 1, seed = 3)
cat("  data dim:", dim(sim_uni$data), "\n")
cat("  Phi:", round(sim_uni$truth$Phi, 3), "\n\n")

cat("=== All tests passed! ===\n")
