#include license.stan

functions {
#include chunks/functions-utils.stan
#include chunks/functions-prior.stan
#include chunks/functions-approx.stan
}

data {
#include chunks/data-common.stan
#include chunks/data-approx.stan
  int y_int[1, num_obs]; // response variable (binary)
}

transformed data{
#include chunks/tdata-approx.stan
  vector[num_obs] glm_a = rep_vector(0.0, num_obs);
}

parameters {
#include chunks/parameters-common.stan
  vector[sum(num_xi)] xi;
}

model {
  vector[sum(num_xi)] glm_b = STAN_build_glm_b(components,
    num_xi, C_ranks, seq_B, L, alpha, ell, xi);
  xi ~ normal(0, 1);
#include chunks/model-priors_common.stan
  // target += normal_lpdf(y_norm | MU, sigma);
  y_int[1] ~ bernoulli_logit_glm(PSI_mats, glm_a, glm_b);
}

//generated quantities {
//  vector[num_obs] f_latent[num_comps] = STAN_build_f(components,
//    num_xi, C_ranks, seq_B, L, PSI_mats, alpha, ell, xi);
//}
