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
functions {
#include functions-utils.stan
#include functions-prior.stan
#include functions-approx.stan
}

data {
#include data.stan
  real<lower=0> X_sd[num_cov_cont];
  real<lower=0> X_hr[num_cov_cont];
  real<lower=0> scale_bf; // factor c to determine domain width L
  int<lower=0> num_bf;    // number of basis functions
  int<lower=0> num_xi[num_comps];
  
  // Categorical decomposition stuff
  int<lower=0> C_sizes[num_comps];
  int<lower=0> C_ranks[num_comps];
  int<lower=0> C_rsp[num_comps]; // element wise product of C_sizes and C_ranks
  vector[sum(C_rsp)] C_vecs;
  real C_vals[sum(C_ranks)];
}

transformed data{
  vector[num_bf] seq_M = STAN_seq_len(num_bf);
  real L[num_cov_cont];
  for(l in 1:num_cov_cont) {
    L[l] = X_hr[l] * scale_bf;
  }
  //matrix[num_obs, num_bf] PHI_mats[num_cov_cont] = 
  //  STAN_create_phi_mats(num_obs, num_cov_cont, seq_M, scale_bf, x_cont, X_hr);
  //matrix[num_obs, sum(num_xi)] PSI_mats = STAN_create_psi_mats(num_obs, num_bf,
  //  num_xi, PHI_mats, components, x_cat, C_vals, C_vecs, C_ranks, C_sizes, C_rsp);
}

parameters {
#include parameters.stan
  vector[num_bf] xi; // sum(num_xi)
}

transformed parameters {
 vector[num_obs] f_latent[1];
 {
   vector[num_bf] dj = STAN_basisfun_eq_multipliers(alpha[1], ell[1], seq_M, L[1]);
   f_latent[1] = STAN_basisfun_eq(x_cont[1], seq_M, L[1]) * (dj .* xi);
 }
      
 //vector[num_obs] f_latent[num_comps] = STAN_build_f_latent(num_obs, components,
 // num_xi, C_ranks, seq_M, X_hr, scale_bf, PSI_mats, alpha, ell, xi);
}

model {
  vector[num_obs] f_sum = STAN_vectorsum(f_latent, num_obs) + c_hat;
  xi ~ normal(0, 1);
#include model.stan
}
