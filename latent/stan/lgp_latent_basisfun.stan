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
  
  // PHI
  real L = 1.0 * scale_bf;
  for(ix in 1:num_cov_cont) {
    for(m in 1:num_bf) {
      PHI_mats[ix][:,m] = STAN_phi(x_cont[ix], m, L);
    }
  }
  
  // VAR_PHI 2
  for(n in 1:num_obs) {
    int z2 = x_cat[idx_z2, n];
    for(c in 1:C2_num) {
      VP2[n,c] = sqrt(C2_vals[c]) * C2_vecs[z2,c];
    }
  }
  
  // VAR_PHI 3
  for(n in 1:num_obs) {
    int z3 = x_cat[idx_z3, n];
    for(c in 1:C3_num) {
      VP3[n,c] = sqrt(C3_vals[c]) * C3_vecs[z3,c];
    }
  }
}

parameters {
  real<lower=1e-12> alpha[num_comps]; // component magnitudes
  real<lower=1e-12> ell[num_ell]; // lengthscales
  real<lower=1e-12> sigma[obs_model==1];
  real<lower=1e-12> phi[obs_model==3];
  real<lower=1e-12, upper=1-1e-12> gamma[obs_model==5];
  
  vector[num_bf] xi_1[1];
  vector[num_bf] xi_2[C2_num];
  vector[num_bf] xi_3[C3_num];
}

transformed parameters {
  vector[num_obs] f_latent[num_comps];
  {
    // Compute diagonals of lambda matrices
    vector[num_bf] d1 = STAN_lambda_matrix(ell[1], num_bf, L);
    vector[num_bf] d2 = STAN_lambda_matrix(ell[2], num_bf, L);
    vector[num_bf] d3 = STAN_lambda_matrix(ell[3], num_bf, L);

    // Build the three components
    f_latent[1] = PHI_mats[1] * (d1 .* xi_1[1]);   // (N,M) x (M)
    f_latent[2] = rep_vector(0, num_obs);
    f_latent[3] = rep_vector(0, num_obs);
    for(m in 1:num_bf){
      vector[num_obs] phi_m = PHI_mats[1][:,m];
      real d_m = sqrt(d2[m]);
      for(c in 1:C2_num) {
        vector[num_obs] v_c = VP2[:,c];
        f_latent[2] += d_m * xi_2[c][m] * phi_m .* v_c;
      }
    }
    for(m in 1:num_bf){
      vector[num_obs] phi_m = PHI_mats[1][:,m];
      real d_m = sqrt(d3[m]);
      for(c in 1:C3_num) {
        vector[num_obs] v_c = VP3[:,c];
        f_latent[3] += d_m * xi_3[c][m] * phi_m .* v_c;
      }
    }
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
