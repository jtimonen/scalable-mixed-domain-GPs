/*
  Author: Juho Timonen
    
  This is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  This code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this code.  If not, see <http://www.gnu.org/licenses/>.
*/

functions{
#include functions-utils.stan
#include functions-prior.stan
#include functions-basisfun.stan
}

data {

  // Dimensions
  int<lower=0> num_obs;         // number of observations
  int<lower=0> num_cov_cont;    // number of continuous covariates
  int<lower=0> num_cov_cat;     // number of categorical covariates
  int<lower=1> num_comps;       // number of additive components
  int<lower=0> num_ell;         // number of lengthscale parameters
  int<lower=0> components[num_comps, 9];
  int<lower=0> x_cat_num_levels[num_cov_cat];

  // Covariates and response variable
  vector[num_obs] x_cont[num_cov_cont];
  int x_cat[num_cov_cat, num_obs];

  // Priors: types and hyperparameters
  int<lower=0> prior_alpha[num_comps, 2];  // {prior_type, transform}
  int<lower=0> prior_ell[num_ell, 2];      // {prior_type, transform}
  real hyper_alpha[num_comps, 3];
  real hyper_ell[num_ell, 3];
  
  // Inputs specific to latent model
  int<lower=1,upper=5> obs_model; // 1-5: Gaussian, Poisson, NB, Bin, BB
  int<lower=0> y_int[obs_model>1, num_obs]; // response variable (int)
  real y_real[obs_model==1, num_obs]; // response variable (real)
  int<lower=1> y_num_trials[obs_model>3, num_obs]; // for Bin or BB model
  int<lower=0> prior_sigma[obs_model==1, 2];
  int<lower=0> prior_phi[obs_model==3, 2];
  real hyper_sigma[obs_model==1, 3];
  real hyper_phi[obs_model==3, 3];
  real hyper_gamma[obs_model==5, 2];
  vector[num_obs] c_hat; // GP mean vector
  
  // Inputs related to basis function approximation
  real<lower=0> scale_bf; // factor c to determine domain width L
  int<lower=1> num_bf;    // number of basis functions
  
  // Categorical decomposition stuff
  int<lower=1> C2_num;
  int<lower=1> C3_num;
  matrix[C2_num, C2_num] C2_vecs;
  matrix[C3_num, C3_num] C3_vecs;
  vector[C2_num] C2_vals;
  vector[C3_num] C3_vals;
  
}

transformed data{
  int idx_z2 = components[2, 8];
  int idx_z3 = components[3, 8];
  matrix[num_obs, num_bf] PHI_mats[num_cov_cont];
  matrix[num_obs, C2_num] VP2;
  matrix[num_obs, C3_num] VP3;
  matrix[num_obs, num_bf] PSI_1;
  matrix[num_obs, num_bf*C2_num] PSI_2;
  matrix[num_obs, num_bf*C3_num] PSI_3;
  vector[num_bf] seq_C1 = STAN_seq_len_rep_times(num_bf, 1);
  vector[num_bf*C2_num] seq_C2 = STAN_seq_len_rep_times(num_bf, C2_num);
  vector[num_bf*C3_num] seq_C3 = STAN_seq_len_rep_times(num_bf, C3_num);
  
  // PHI
  real L = 1.0 * scale_bf;
  for(ix in 1:num_cov_cont) {
    PHI_mats[ix] = STAN_PHI_eq(x_cont[ix], seq_C1, L);
  }
  
  // PSI_1
  PSI_1 = PHI_mats[1];
  
  // PSI_2
  for(c in 1:C2_num) {
    int i1 = 1 + (c-1)*num_bf;
    int i2 = c*num_bf;
    vector[num_obs] vpc = sqrt(C2_vals[c]) * C2_vecs[x_cat[idx_z2, :],c];
    PSI_2[:, i1:i2] = PHI_mats[1] .* rep_matrix(vpc, num_bf);
  }
  
    // PSI_3
  for(c in 1:C3_num) {
    int i1 = 1 + (c-1)*num_bf;
    int i2 = c*num_bf;
    vector[num_obs] vpc = sqrt(C3_vals[c]) * C3_vecs[x_cat[idx_z3, :],c];
    PSI_3[:, i1:i2] = PHI_mats[1] .* rep_matrix(vpc, num_bf);
  }
  
}

parameters {
  real<lower=1e-12> alpha[num_comps]; // component magnitudes
  real<lower=1e-12> ell[num_ell]; // lengthscales
  real<lower=1e-12> sigma[obs_model==1];
  real<lower=1e-12> phi[obs_model==3];
  real<lower=1e-12, upper=1-1e-12> gamma[obs_model==5];
  
  vector[num_bf] xi_1;
  vector[num_bf*C2_num] xi_2;
  vector[num_bf*C3_num] xi_3;
}

transformed parameters {
  vector[num_obs] f_latent[num_comps];
  {
    
    // Multipliers
    vector[num_bf] d1 = STAN_diag_spd_eq(alpha[1], ell[1], seq_C1, L);
    vector[num_bf*C2_num] d2 = STAN_diag_spd_eq(alpha[2], ell[2], seq_C2, L);
    vector[num_bf*C3_num] d3 = STAN_diag_spd_eq(alpha[3], ell[3], seq_C3, L);
    
    // Build the three components
    f_latent[1] = PSI_1 * (d1 .* xi_1);
    f_latent[2] = PSI_2 * (d2 .* xi_2);
    f_latent[3] = PSI_3 * (d3 .* xi_3);
  }
}

model {
  vector[num_obs] f_sum = STAN_vectorsum(f_latent, num_obs) + c_hat;
  
  // Multiplier priors
  xi_1[1] ~ normal(0, 1);
  for(c in 1:C2_num) {
    xi_2[c] ~ normal(0, 1);
  }
  for(c in 1:C3_num) {
    xi_3[c] ~ normal(0, 1);
  }
  
  // Parameter priors
  for(j in 1:num_comps){
    target += STAN_log_prior(alpha[j], prior_alpha[j], hyper_alpha[j]);
  }
  for(j in 1:num_ell){
    target += STAN_log_prior(ell[j], prior_ell[j], hyper_ell[j]);
  }
  if(obs_model==1){
    target += STAN_log_prior(sigma[1], prior_sigma[1], hyper_sigma[1]);
  }else if(obs_model==3){
    target += STAN_log_prior(phi[1], prior_phi[1], hyper_phi[1]);
  }else if(obs_model==5){
    target += beta_lpdf(gamma[1] | hyper_gamma[1][2], hyper_gamma[1][2]);
  }

  // Likelihood
  if(obs_model==1) {
    // 1. Gaussian
    real MU[num_obs] = to_array_1d(f_sum); // means
    target += normal_lpdf(y_real[1] | MU, sigma[1]);
  }else if(obs_model==2){
    // 2. Poisson
    real LOG_MU[num_obs] = to_array_1d(f_sum); // means (log-scale)
    target += poisson_log_lpmf(y_int[1] | LOG_MU);
  }else if(obs_model==3){
    // 3. Negative binomial
    real LOG_MU[num_obs] = to_array_1d(f_sum); // means (log-scale)
    real PHI[num_obs] = to_array_1d(rep_vector(phi[1], num_obs)); // dispersion
    target += neg_binomial_2_log_lpmf(y_int[1] | LOG_MU, PHI);
  }else if(obs_model==4){
    // 4. Binomial
    real LOGIT_P[num_obs] = to_array_1d(f_sum); // p success (logit-scale)
    target += binomial_logit_lpmf(y_int[1] | y_num_trials[1], LOGIT_P);
  }else if(obs_model==5){
    // 5. Beta-binomial
    real tgam = inv(gamma[1]) - 1.0;
    vector[num_obs] P = inv_logit(f_sum); // p success
    real aa[num_obs] = to_array_1d(P * tgam);
    real bb[num_obs] = to_array_1d((1.0 - P) * tgam);
    target += beta_binomial_lpmf(y_int[1] | y_num_trials[1], aa, bb);
  }
}
