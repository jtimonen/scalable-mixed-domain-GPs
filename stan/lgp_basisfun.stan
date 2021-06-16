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
#include functions-prior.stan
#include functions-kernels.stan
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
  real delta; // jitter to ensure pos. def. kernel matrices

  // Categorical and continuous covariates, and target variable y
  vector[num_obs] x_cont[num_cov_cont];
  int x_cat[num_cov_cat, num_obs];
  vector[num_obs] y_norm;

  // Priors: types and hyperparameters
  int<lower=0> prior_alpha[num_comps, 2];  // {prior_type, transform}
  int<lower=0> prior_ell[num_ell, 2];      // {prior_type, transform}
  real hyper_alpha[num_comps, 3];
  real hyper_ell[num_ell, 3];
  int<lower=0> prior_sigma[1, 2];
  real hyper_sigma[1, 3];
  
  // Inputs related to basis function approximation
  real<lower=0> c_bf;   // factor c to determine domain width L
  int<lower=1> M_bf;    // number of basis functions
}

transformed data{
  
  // Constant kernel matrix ranks and eigenvalues
  int ranks[num_comps] = STAN_ranks(components, x_cat_num_levels);
  int R = sum(ranks);
  int RM = R * num_basisfun;

  // Compute L and PHI for each component
  real L_bf[num_comps];
  matrix[num_obs, M_bf] PHI[num_comps];
  for(j in 1:num_comps) {
    int idx_cont = components[j, 9];
    print("idx_cont=", idx_cont);
    //L_bf[j] = c_bf*max(x_cont[idx_cont]);
    //PHI[j] = PHI_EQ(num_obs, M_bf, L_bf[j], x_cont[idx_cont]);
  }
  
}

parameters {
  real<lower=1e-12> alpha[num_comps]; // component magnitudes
  real<lower=1e-12> ell[num_ell]; // lengthscales
  real<lower=1e-12> sigma[1]; // noise std
  vector[M_bf] xi[num_components]; # basis function coefficients
}


transformed parameters {
  // Compute spectral densities
  vector[num_obs] f[num_comps];
  for (j in 1:num_comps){
    vector[M_bf] diagSPD_j = diagSPD_EQ(alpha[j], ell[j], L_bf[j], M_);
    f[j] = PHI_[j] * (diagSPD_j .* beta_f)
  }
}

model {
      
  // Likelihood
  vector[num_obs] f_sum = rep_vector(0.0, num_obs);
  for(j in 1:num_comps){
    f_sum += f[j];
  }
  y_norm ~ normal(f_sum, sigma[1]);
  
  // Priors
  for(j in 1:num_comps){
    target += STAN_log_prior(alpha[j], prior_alpha[j], hyper_alpha[j]);
  }
  for(j in 1:num_ell){
    target += STAN_log_prior(ell[j], prior_ell[j], hyper_ell[j]);
  }
  target += STAN_log_prior(sigma[1], prior_sigma[1], hyper_sigma[1]);
}

