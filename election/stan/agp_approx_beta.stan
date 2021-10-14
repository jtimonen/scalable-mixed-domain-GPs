#include license.stan

functions {
#include chunks/functions-utils.stan
#include chunks/functions-prior.stan
#include chunks/functions-approx.stan
}

data {
#include chunks/data-common.stan
#include chunks/data-approx.stan
  int<lower=0> votes_rep[num_obs];
  int<lower=0> votes_dem[num_obs];
  int<lower=1> votes_total[num_obs];
  vector<lower=0,upper=1>[num_obs] rep_share;
}

transformed data{
#include chunks/tdata-approx.stan
  vector[num_obs] glm_a = rep_vector(0.0, num_obs);
}

parameters {
#include chunks/parameters-common.stan
  vector[sum(num_xi)] xi;
  real<lower=0.000000001> gamma;
  real mu;
}

model {
  vector[sum(num_xi)] glm_b = STAN_build_glm_b(components,
    num_xi, C_ranks, seq_B, L, alpha, ell, xi);
  //vector[num_obs] obs_mu = nu*inv_logit(mu + PSI_mats *  glm_b);
  xi ~ normal(0, 1);
  gamma ~ lognormal(1.0, 1.0);
  mu ~ normal(0, 0.5);
#include chunks/model-priors_common.stan
  //rep_share ~ beta(obs_mu, nu-obs_mu);
  {
    real tgam = inv(gamma) - 1.0;
    vector[num_obs] P = inv_logit(mu + PSI_mats *  glm_b); // p success
    vector[num_obs] aa = P * tgam;
    vector[num_obs] bb = (1.0 - P) * tgam;
    target += beta_binomial_lpmf(votes_rep | votes_total, aa, bb);
  }

}

//generated quantities {
//  vector[num_obs] f_latent[num_comps] = STAN_build_f(components,
//    num_xi, C_ranks, seq_B, L, PSI_mats, alpha, ell, xi);
//}
