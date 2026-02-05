sim_mlVAR <- function(
    N,
    T_obs,
    p,
    q,
    seed) {
  # covariates
  for (variable in vector) {

  }
}


set.seed(1337)

N     <- 25
T_obs <- 100
p     <- 4      # outcomes y_1..y_4
q     <- 2      # covariates x_1, x_2

# ----------------------------
# 1) Covariates (person-time)
# ----------------------------
# x_1: time-varying continuous, x_2: binary or continuous (choose)
x_1 <- matrix(rnorm(N * T_obs, 0, 1), N, T_obs)
x_2 <- matrix(rbinom(N * T_obs, 1, 0.5), N, T_obs)  # binary covariate
X   <- array(NA_real_, dim = c(N, T_obs, q))
X[,,1] <- x_1
X[,,2] <- x_2

# ----------------------------
# 2) Population-level params
# ----------------------------
alpha_pop <- c(-0.5, -0.2, -0.7, -0.4)         # baseline log-odds per node
gamma_pop <- matrix(c( 0.6, -0.4,              # effects of x_1, x_2 on y1
                       0.3,  0.2,              # on y2
                       -0.2,  0.5,              # on y3
                       0.1, -0.3),             # on y4
                    nrow = p, byrow = TRUE)

# VAR(1) population-level edge matrix (rows: target node j, cols: lagged k)
B_pop <- matrix(0, p, p)
diag(B_pop) <- c(1.0, 0.8, 0.9, 0.7)           # autoregressive self-effects

# some cross-lag effects (keep modest for stability)
B_pop[2,1] <-  0.4
B_pop[1,2] <- -0.3
B_pop[3,2] <-  0.2
B_pop[4,3] <- -0.25

# ----------------------------
# 3) Random-effects scales
# ----------------------------
sd_alpha <- 0.6     # between-person SD in intercepts
sd_B     <- 0.25    # between-person SD per edge (iid here for simplicity)

# Draw person-specific intercepts and VAR matrices
alpha_i <- matrix(NA_real_, N, p)
B_i     <- array(NA_real_, dim = c(N, p, p))

for (i in 1:N) {
  alpha_i[i,] <- rnorm(p, mean = alpha_pop, sd = sd_alpha)
  B_i[i,,]    <- B_pop + matrix(rnorm(p*p, 0, sd_B), p, p)
}

inv_logit <- function(z) 1/(1 + exp(-z))

# ----------------------------
# 4) Simulate binary outcomes
# ----------------------------
Y <- array(0L, dim = c(N, T_obs, p))

# Initialize t=1 (could also set all zeros or draw from intercept-only probs)
for (i in 1:N) {
  eta1 <- alpha_i[i,] + X[i,1,] %*% t(gamma_pop)  # 1 x p
  p1   <- inv_logit(as.numeric(eta1))
  Y[i,1,] <- rbinom(p, 1, p1)
}

# Forward simulate
for (i in 1:N) {
  for (t in 2:T_obs) {
    lag_y <- Y[i,t-1,]                      # length p
    # linear predictor for each node j:
    eta <- alpha_i[i,] +
      as.numeric(X[i,t,] %*% t(gamma_pop)) + # covariates -> p-vector
      as.numeric(B_i[i,,] %*% lag_y)         # VAR(1) -> p-vector

    prob <- inv_logit(eta)
    Y[i,t,] <- rbinom(p, 1, prob)
  }
}

# ----------------------------
# 5) Put into a long data frame
# ----------------------------
df <- do.call(rbind, lapply(1:N, function(i) {
  data.frame(
    id = i,
    t  = 1:T_obs,
    x_1 = X[i,,1],
    x_2 = X[i,,2],
    y_1 = Y[i,,1],
    y_2 = Y[i,,2],
    y_3 = Y[i,,3],
    y_4 = Y[i,,4]
  )
}))
