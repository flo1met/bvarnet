// TODO:
// - missing data skipping of lags (handle inside missing data handling R function) - DONE
// - different priors - DONE
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

    // §1.5: partial sum for reduce_sum parallelisation
    real partial_sum_binary(
        array[] int y_slice,
        int start, int end,
        matrix X_fixed,
        vector eta_re,
        vector b_fixed
    ) {
        return bernoulli_logit_glm_lpmf(y_slice |
            X_fixed[start:end], eta_re[start:end], b_fixed);
    }
}

data {
    int<lower=1> p; // nb outcome parameters
    //int<lower=0> w; // nb interactions

    int<lower=1> J; // nb persons/groups
    int<lower=1> K; // AR(K)

    int<lower=1> n_obs; // legth Y, J*T-K (-K bc 1 row will be deleted as no lagged effects, therefore no estimation, when we add indicator variable this will be removed(?))
    int<lower=1> n_fe; // number of fixed effects, columns of X, max q+w+1
    int<lower=0> n_re; // number of random effects, columns of Z, max q+p+w+1

    array[n_obs] int<lower=1, upper=J> id; // person identifier

    array[n_obs, p] int<lower=0, upper=1> Y; // outcome matrix in long format

    matrix[n_obs, p*K] B; // design matrix lagged parameters
    matrix[n_obs, n_fe] X; // design matrix fixed effects (centered!)
    matrix[n_obs, n_re] Z; // design matrix random effects

    /// PRIORS: 1=normal, 2=t, 3=chauchy
    // beta
    int<lower=1> prior_beta_fam; 
    real beta_loc;
    real<lower=0> beta_scale;
    real<lower=0> beta_df;

    int<lower=1> prior_phi_fam;
    real phi_loc;
    real<lower=0> phi_scale;
    real<lower=0> phi_df;

    int<lower=1> prior_sd_fam;
    real sd_loc;
    real<lower=0> sd_scale;
    real<lower=0> sd_df;

    int<lower=1> grainsize; // reduce_sum grain size (1 = auto)
}
transformed data {
   matrix[n_obs, n_fe + p*K] X_fixed = append_col(X,B); // intercept correct? yes, as long as its in X!
}
parameters {
   matrix[p*K, p] phi; // lag coefficient matrix. Flat for glm 
   matrix[n_fe, p] beta; // fixed effects

   array[p] matrix[J, n_re] z_u; // non-centered parametrisation, latent parameter
   matrix<lower=0>[p, n_re] sd_u; // RE scales
}
transformed parameters {
   // random effects precomputation // optimise regarding loops and vectorisation
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

    /// std_normal prior for latent mean
    for (node in 1:p)
        target += std_normal_lpdf(to_vector(z_u[node]));

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

        // §1.5: parallelise bernoulli_logit_glm over observations
        target += reduce_sum(partial_sum_binary, Y[, node],
                             grainsize, X_fixed, eta_re, b_fixed);
    }

}
