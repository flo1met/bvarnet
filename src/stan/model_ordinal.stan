// Adjacent-category ordinal model
// Reduces to binary logistic regression when C = 2
// See dev/model_development_plan.md §2 for specification & rationale
//
// TODO:
// - per-node varying C (ragged cutpoints)
// - dummy-coded ordinal lags
// - reduce_sum parallelisation
// - correlated random effects

#include /functions.stan

data {
    int<lower=1> p; // nb outcome parameters

    int<lower=1> J; // nb persons/groups
    int<lower=1> K; // AR(K)

    int<lower=1> n_obs; // number of modeled observations
    int<lower=1> n_fe; // number of fixed effects (columns of X, NO intercept for ordinal)
    int<lower=0> n_re; // number of random effects (columns of Z)

    array[n_obs] int<lower=1, upper=J> id; // person identifier

    // --- ordinal-specific ---
    int<lower=2> C; // number of ordered categories (shared across all nodes) // keeping fixed for now
    array[n_obs, p] int<lower=1, upper=C> Y; // ordinal outcome matrix

    matrix[n_obs, p*K] B; // design matrix lagged parameters
    matrix[n_obs, n_fe] X; // design matrix fixed effects (NO intercept — absorbed by kappa)
    matrix[n_obs, n_re] Z; // design matrix random effects

    /// PRIORS: 1=normal, 2=t, 3=cauchy
    // beta
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

    // kappa (cutpoints)
    int<lower=1> prior_kappa_fam;
    real kappa_loc;
    real<lower=0> kappa_scale;
    real<lower=0> kappa_df;
}
transformed data {
   matrix[n_obs, n_fe + p*K] X_fixed = append_col(X, B);
}
parameters {
   matrix[p*K, p] phi; // lag coefficient matrix
   matrix[n_fe, p] beta; // fixed effects (no intercept row)

   array[p] matrix[J, n_re] z_u; // non-centered parametrisation, latent parameter
   matrix<lower=0>[p, n_re] sd_u; // RE scales

   // --- ordinal-specific ---
   array[p] ordered[C - 1] kappa; // cutpoints per node (ordered enforces kappa_1 < kappa_2 < ...)
}
transformed parameters {
   // random effects precomputation (non-centered)
   array[p] matrix[J, n_re] u;
   for (node in 1:p)
    for (j in 1:J)
        for (k in 1:n_re)
            u[node][j,k] = sd_u[node,k] * z_u[node][j,k];
}
model {
    /// priors
    target += set_prior(to_vector(beta), prior_beta_fam, beta_loc, beta_scale, beta_df);
    target += set_prior(to_vector(phi), prior_phi_fam, phi_loc, phi_scale, phi_df);

    if (n_re > 0)
        target += set_half_prior(to_vector(sd_u), prior_sd_fam, sd_loc, sd_scale, sd_df);

    /// std_normal prior for latent mean (non-centered RE)
    for (node in 1:p)
        target += std_normal_lpdf(to_vector(z_u[node]));

    /// cutpoint priors
    for (node in 1:p)
        target += set_prior(kappa[node], prior_kappa_fam, kappa_loc, kappa_scale, kappa_df);

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

        // linear predictor (no intercept — absorbed into kappa)
        vector[n_obs] eta = X_fixed * b_fixed + eta_re;

        // adjacent-category likelihood via categorical_logit_lpmf
        // Pre-compute cumulative sum of kappa (constant across observations)
        vector[C] kappa_cumsum;
        kappa_cumsum[1] = 0;
        for (c in 2:C)
            kappa_cumsum[c] = kappa_cumsum[c - 1] + kappa[node][c - 1];

        for (i in 1:n_obs) {
            vector[C] lambda;
            for (c in 1:C)
                lambda[c] = (c - 1) * eta[i] - kappa_cumsum[c];
            target += categorical_logit_lpmf(Y[i, node] | lambda);
        }
    }
}
