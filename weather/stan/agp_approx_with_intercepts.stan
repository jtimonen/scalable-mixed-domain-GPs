#include license.stan

functions {
#include chunks/functions-utils.stan
#include chunks/functions-prior.stan
#include chunks/functions-approx.stan
}

data {
#include chunks/data-common.stan
#include chunks/data-approx.stan
  real y_norm[num_obs]; // normalized response variable (real)
  int<lower=0> prior_sigma[1, 2];
  real hyper_sigma[1, 3];
  int<lower=0> num_ids;
  int<lower=1, upper=num_ids> idx_expand_intercepts[num_obs];
}

transformed data{
#include chunks/tdata-approx.stan
}

parameters {
#include chunks/parameters-common.stan
  real<lower=1e-12> sigma;
  vector[sum(num_xi)] xi;
  vector[num_ids] id_intercepts;
}

transformed parameters {
  vector[num_obs] f_latent[num_comps] = STAN_build_f(components,
    num_xi, C_ranks, seq_B, L, PSI_mats, alpha, ell, xi);
}

model {
  real MU[num_obs] = to_array_1d(STAN_vectorsum(f_latent, num_obs) +
      id_intercepts[idx_expand_intercepts]);
  xi ~ normal(0, 1);
  target += STAN_log_prior(sigma, prior_sigma[1], hyper_sigma[1]);
#include chunks/model-priors_common.stan
  id_intercepts ~ normal(0.0, 1.0); // prior for intercepts
  target += normal_lpdf(y_norm | MU, sigma);
}
