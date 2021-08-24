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
  int<lower=1> len_eigvals;
  int<lower=1> len_eigvecs;
  int<lower=0> C_sizes[num_comps];
  vector<lower=0>[len_eigvals] C_eigvals;
  vector[len_eigvecs] C_eigvecs;
  int<lower=0> C_inds[num_comps, 2];
  
}

transformed data{
  vector[N] PHI[num_cov_cont, num_bf];
  matrix[N] VAR_PHI[len_eigvals];
  int idx_cc = 0;
  
  // PHI
  real L = 1.0 * scale_bf;
  for(ix in 1:num_cov_cont) {
    for(m in 1:num_bf) {
      PHI[ix,m] = STAN_phi(x_cont[ix], m, L);
    }
  }
  
  // VAR_PHI
  for(j in 1:num_comps) {
    int i1;
    int i2;
    if(components[j,1] != 1) {
      idx_cc += 1;
      int i1 = C_inds[j, 1];
      int i2 = C_inds[j, 2];
    }
    for(m in 1:num_bf) {
      PHI[ix,m] = STAN_phi(x_cont[ix], m, L);
    }
  }
}

parameters {
  real<lower=1e-12> alpha[num_comps]; // component magnitudes
  real<lower=1e-12> ell[num_ell]; // lengthscales
  real<lower=1e-12> sigma[obs_model==1];
  real<lower=1e-12> phi[obs_model==3];
  real<lower=1e-12, upper=1-1e-12> gamma[obs_model==5];
  
  vector[num_bf] xi[num_comps, len_eigvals]; // basis function multipliers
}

transformed parameters {
  vector[num_obs] f_latent[num_comps];
  {
    matrix[N, num_bf] PSI[num_comps];
    int idx_ell = 0;
    int idx_alpha = 0;
  
    // Loop through components
    for(j in 1:num_comps){
      

      
      // 1. Initialize with constant part of PSI
      matrix[N, num_bf] PSI_j = rep_matrix(1.0, N, num_bf);
      
      // 2. Get component properties
      int opts[9] = components[j];
      int ctype = opts[1]; # changes in 2.0
      int idx_cont = opts[9]; # changes in 2.0
      
      // 3. Pick the possible continuous covariate of this component
      if (ctype > 0) {
        PSI_j = PSI_j .* PHI[idx_cont];
      }
      
      // 4. Pick the possible categorical covariate of this component
      if (ctype != 1) {
        PSI_j = PSI_j .* PHI[idx_cont];
      }
      
      KX[j] = K; // store kernel matrix
    }
    return(KX);
    
    for(j in 1:num_comps){
      f_latent[j] = intercept + PHI_f * (diagSPD_f .* beta_f);
    }
  }
}

model {
  xi ~ normal(0, 1);

  vector[num_obs] f_sum = STAN_vectorsum(f_latent, num_obs) + c_hat;
  for(j in 1:num_comps){ target += normal_lpdf(eta[j] | 0, 1);}
  
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
