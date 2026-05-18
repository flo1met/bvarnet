// Gaussian (continuous, identity link) model
// See dev/model_development_plan.md §3 for specification & rationale
//
// TODO:
// - reduce_sum parallelisation
// - correlated random effects

functions {
    real set_prior(vector x, int fam, real loc, real scale, real df) {
        if (fam == 1) return normal_lpdf(x | loc, scale);
        if (fam == 2) return student_t_lpdf(x | df, loc, scale);
        if (fam == 3) return cauchy_lpdf(x | loc, scale);
        reject("Unknown prior family code:", fam);
    }

    real set_half_prior(vector x, int fam, real loc, real scale, real df) {
        if (fam == 1) return normal_lpdf(x | loc, scale) - rows(x) * normal_lccdf(0 | loc, scale);
        if (fam == 2) return student_t_lpdf(x | df, loc, scale) - rows(x) * student_t_lccdf(0 | df, loc, scale);
        if (fam == 3) return cauchy_lpdf(x | loc, scale) - rows(x) * cauchy_lccdf(0 | loc, scale);
        reject("Unknown prior family code:", fam);
   }
}

data {
    int<lower=1> p; // nb outcome parameters

    int<lower=1> J; // nb persons/groups
    int<lower=1> K; // AR(K)

    int<lower=1> n_obs; // number of modeled observations
    int<lower=1> n_fe; // number of fixed effects (columns of X, incl intercept)
    int<lower=0> n_re; // number of random effects (columns of Z)

    array[n_obs] int<lower=1, upper=J> id; // person identifier

    // --- gaussian-specific ---
    matrix[n_obs, p] Y; // real-valued outcome matrix

    matrix[n_obs, p*K] B; // design matrix lagged parameters
    matrix[n_obs, n_fe] X; // design matrix fixed effects (with intercept)
    matrix[n_obs, n_re] Z; // design matrix random effects

    /// PRIORS: 1=normal, 2=t, 3=cauchy
    // intercept (beta row 1)
    int<lower=1> prior_intercept_fam;
    real intercept_loc;
    real<lower=0> intercept_scale;
    real<lower=0> intercept_df;

    // beta (slopes: beta rows 2..n_fe)
    int<lower=1> prior_beta_fam;
    real beta_loc;
    real<lower=0> beta_scale;
    real<lower=0> beta_df;

    // phi
    int<lower=1> prior_phi_fam;
    real phi_loc;
    real<lower=0> phi_scale;
    real<lower=0> phi_df;

    // sd_u
    int<lower=1> prior_sd_fam;
    real sd_loc;
    real<lower=0> sd_scale;
    real<lower=0> sd_df;

    // sigma (residual SD)
    int<lower=1> prior_sigma_fam;
    real sigma_loc;
    real<lower=0> sigma_scale;
    real<lower=0> sigma_df;
}
transformed data {
   matrix[n_obs, n_fe + p*K] X_fixed = append_col(X, B);
}
parameters {
   matrix[p*K, p] phi; // lag coefficient matrix
   matrix[n_fe, p] beta; // fixed effects (incl intercept row)

   array[p] matrix[J, n_re] z_u; // non-centered parametrisation, latent parameter
   matrix<lower=0>[p, n_re] sd_u; // RE scales

   // --- gaussian-specific ---
   vector<lower=0>[p] sigma; // per-node residual SD
}
transformed parameters {
   // random effects precomputation (non-centered)
   array[p] matrix[J, n_re] u;
   for (node in 1:p)
        u[node] = z_u[node] .* rep_matrix(sd_u[node,], J);
}
model {
    /// priors
    target += set_prior(to_vector(beta[1,]), prior_intercept_fam, intercept_loc, intercept_scale, intercept_df);
    if (n_fe > 1)
        target += set_prior(to_vector(beta[2:n_fe,]), prior_beta_fam, beta_loc, beta_scale, beta_df);
    target += set_prior(to_vector(phi), prior_phi_fam, phi_loc, phi_scale, phi_df);

    if (n_re > 0)
        target += set_half_prior(to_vector(sd_u), prior_sd_fam, sd_loc, sd_scale, sd_df);

    /// std_normal prior for latent mean (non-centered RE)
    for (node in 1:p)
        target += std_normal_lpdf(to_vector(z_u[node]));

    /// residual SD prior (half-prior)
    target += set_half_prior(sigma, prior_sigma_fam, sigma_loc, sigma_scale, sigma_df);

    // Likelihood
    for (node in 1:p) {
        vector[n_obs] eta_re;
        vector[n_fe + p*K] b_fixed;

        // calculate offset: eta_re[n] = Z[n,] * u[node][id[n], ]'
        eta_re = (n_re > 0)
            ? rows_dot_product(Z, u[node][id,])
            : rep_vector(0.0, n_obs);

        // get combined FE vector
        b_fixed[1:n_fe] = beta[,node];
        b_fixed[(n_fe + 1):(n_fe + p*K)] = phi[,node];

        // GLM form: normal_id_glm fuses X*b + offset into one kernel
        target += normal_id_glm_lpdf(Y[,node] | X_fixed, eta_re, b_fixed, sigma[node]);
    }
}
