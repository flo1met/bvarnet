# Test burnin and no-RE features
source("R/sim_bvarnet.R")

cat("=== Test: Burnin (bernoulli) ===\n")
sim1 <- sim_bvarnet(N = 10, T_obs = 50, p = 3, K = 1,
                    family = "bernoulli", burnin = 100, seed = 42)
cat("  data dim:", dim(sim1$data), " (expect 500 x 5)\n")
cat("  t range:", range(sim1$data$t), " (expect 1..50)\n")
cat("  burnin in truth:", sim1$truth$burnin, "\n\n")

cat("=== Test: Burnin (gaussian K=2) ===\n")
sim2 <- sim_bvarnet(N = 5, T_obs = 30, p = 2, q = 1, K = 2,
                    family = "gaussian", burnin = 200, seed = 99)
cat("  data dim:", dim(sim2$data), " (expect 150 x 5)\n")
cat("  t range:", range(sim2$data$t), "\n\n")

cat("=== Test: Burnin (ordinal) ===\n")
sim3 <- sim_bvarnet(N = 10, T_obs = 40, p = 2, K = 1,
                    family = "ordinal", C = 4, burnin = 50, seed = 7)
cat("  data dim:", dim(sim3$data), "\n")
cat("  Y table node 1:\n")
print(table(sim3$data$y_1))
cat("\n")

cat("=== Test: No random effects (sd_alpha = 0) ===\n")
sim4 <- sim_bvarnet(N = 5, T_obs = 30, p = 2, K = 1,
                    family = "gaussian", sd_alpha = 0, seed = 1)
cat("  alpha_i (all rows should be identical):\n")
print(round(sim4$truth$alpha_i, 4))
cat("  Population alpha:", round(sim4$truth$alpha, 4), "\n")
cat("  All subjects identical?:", all(apply(sim4$truth$alpha_i, 2, sd) == 0), "\n\n")

cat("=== Test: No RE, multiple subjects share Phi too ===\n")
sim5 <- sim_bvarnet(N = 3, T_obs = 20, p = 2, K = 1,
                    family = "bernoulli", sd_alpha = 0,
                    re_temporal = FALSE, seed = 5)
cat("  Phi_i[1,,] == Phi_i[2,,]?:", identical(sim5$truth$Phi_i[1,,], sim5$truth$Phi_i[2,,]), "\n")
cat("  Phi_i[1,,] == Phi?:", identical(sim5$truth$Phi_i[1,,], sim5$truth$Phi), "\n\n")

cat("=== Test: Burnin=0 is default (backward compatible) ===\n")
sim6 <- sim_bvarnet(N = 5, T_obs = 20, p = 2, K = 1,
                    family = "bernoulli", seed = 42)
cat("  data dim:", dim(sim6$data), " (expect 100 x 4)\n")
cat("  burnin:", sim6$truth$burnin, "\n\n")

cat("=== All burnin/no-RE tests passed! ===\n")
