// Adjacent-category ordinal model
// Reduces to binary logistic regression when C = 2
// See dev/model_development_plan.md §2 for specification & rationale
//
// TODO:
// - per-node varying C (ragged cutpoints)
// - dummy-coded ordinal lags
// - reduce_sum parallelisation (§2.3)
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

   // Fixed 1×C coefficient matrix encoding the adjacent-category constraint:
   // beta_adj[1, c] = c - 1  so that  alpha[c] + eta[n] * (c-1) = lambda[n,c]
   // Passed as the `beta` argument to categorical_logit_glm_lpmf.
   matrix[1, C] beta_adj;
   for (c in 1:C)
       beta_adj[1, c] = c - 1;
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
        u[node] = z_u[node] .* rep_matrix(sd_u[node,], J);
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

        // adjacent-category likelihood via categorical_logit_glm_lpmf
        // (same approach as MaartenMarsman/mixedGM inst/stan/mixed_mrf_conditional.stan)
        //
        // categorical_logit_glm_lpmf(y | X, alpha, beta) computes:
        //   log-prob for obs n, category c  =  alpha[c] + X[n,1] * beta[1,c]
        //
        // With:
        //   X[n, 1]    = eta[n]          (N×1 predictor matrix)
        //   alpha[c]   = -kappa_cumsum[c]
        //   beta[1, c] = c - 1           (fixed in transformed data)
        //
        // => alpha[c] + eta[n]*(c-1) = (c-1)*eta[n] - kappa_cumsum[c] = lambda[n,c]  ✓
        //
        // Pre-compute cumulative kappa
        vector[C] kappa_cumsum;
        kappa_cumsum[1] = 0;
        for (c in 2:C)
            kappa_cumsum[c] = kappa_cumsum[c - 1] + kappa[node][c - 1];

        vector[C] alpha_cat = -kappa_cumsum;
        alpha_cat[1] = 0; // kappa_cumsum[1] = 0 so this is a no-op, but explicit

        target += categorical_logit_glm_lpmf(Y[, node] | to_matrix(eta), alpha_cat, beta_adj);
    }
}
