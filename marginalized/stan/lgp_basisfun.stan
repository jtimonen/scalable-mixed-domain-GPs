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
  real<lower=0> L_bf;    // domain width
  int<lower=1> M_bf;    // number of basis functions
}

transformed data{
  
  // Precompute fixed kernel matrices
  vector[num_obs] m0 = rep_vector(0.0, num_obs);
  matrix[num_obs, num_obs] K_const[num_comps] = STAN_kernel_const_all(
    num_obs, num_obs, x_cat, x_cat, x_cat_num_levels, components
  );
  vector[num_obs] delta_vec = rep_vector(delta, num_obs);

  // Constant kernel matrix ranks and eigenvalues
  int ranks[num_comps] = STAN_ranks(components, x_cat_num_levels);
  int R = sum(ranks);
  int RM = R * M_bf;
  vector[R] bfa_delta = STAN_delta_matrix(K_const, ranks, components);
  matrix[num_obs, R] bfa_theta = STAN_theta_matrix(K_const, ranks, components);
}

parameters {
  real<lower=1e-12> alpha[num_comps]; // component magnitudes
  real<lower=1e-12> ell[num_ell]; // lengthscales
  real<lower=1e-12> sigma[1]; // noise std
}


model {

  // Compute Phi and Lambda
  vector[M_bf] bfa_lambda[num_comps] = STAN_lambda_matrix(ell, 
      M_bf, L_bf, components);
  matrix[num_obs, M_bf] bfa_phi[num_comps] = STAN_phi_matrix(
      x_cont, M_bf, L_bf, components);
  vector[RM] bfa_D = STAN_D_matrix(alpha, bfa_lambda, bfa_delta, ranks); // beta?
  matrix[num_obs, RM] bfa_V = STAN_V_matrix(bfa_phi, bfa_theta, ranks);

  // Approximate likelihood
  target += STAN_multi_normal_bfa_logpdf(y_norm, bfa_V, bfa_D, sigma[1]);
  
  // Priors
  for(j in 1:num_comps){
    target += STAN_log_prior(alpha[j], prior_alpha[j], hyper_alpha[j]);
  }
  for(j in 1:num_ell){
    target += STAN_log_prior(ell[j], prior_ell[j], hyper_ell[j]);
  }
  target += STAN_log_prior(sigma[1], prior_sigma[1], hyper_sigma[1]);
}

