#include license.stan

functions {
#include chunks/functions-utils.stan
#include chunks/functions-prior.stan
#include chunks/functions-approx.stan
}

data {
#include chunks/data-common.stan
#include chunks/data-approx.stan
#include chunks/data-nongaussian.stan
}

transformed data{
#include chunks/tdata-approx.stan
}

parameters {
#include chunks/parameters-common.stan
#include chunks/parameters-nongaussian.stan
  vector[sum(num_xi)] xi;
}

transformed parameters {
 vector[num_obs] f_latent[num_comps] = STAN_build_f(components,
  num_xi, C_ranks, seq_B, L, PSI_mats, alpha, ell, xi);
}

model {
  vector[num_obs] f_sum = STAN_vectorsum(f_latent, num_obs) + c_hat;
  xi ~ normal(0, 1);
#include chunks/model-priors_common.stan
#include chunks/model-nongaussian.stan
}
