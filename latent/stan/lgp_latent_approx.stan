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
  vector<lower=0>[num_cov_cont] X_sd;
  vector<lower=0>[num_cov_cont] X_hr;
  real<lower=0> scale_bf; // factor c to determine domain width L
  int<lower=1> num_bf;    // number of basis functions
  int<lower=0> num_xi[num_comps];
  
  // Categorical decomposition stuff
  int<lower=0> C_sizes[num_comps];
  int<lower=0> C_ranks[num_comps];
  int<lower=0> C_rsp[num_comps]; // element wise product of C_sizes and C_ranks
  vector[sum(C_rsp)] C_vecs;
  real C_vals[sum(C_ranks)];
}

transformed data{
  vector[num_bf] seq_B = STAN_seq_len(num_bf);
  matrix[num_obs, num_bf] mat_B = 
    transpose(rep_matrix(seq_B, num_obs));
  vector[num_cov_cont] L = scale_bf * X_hr;
  matrix[num_obs, num_bf] PHI_mats[num_cov_cont] = 
    STAN_create_basisfun_mats(x_cont, mat_B, L);
  matrix[num_obs, sum(num_xi)] PSI_mats = STAN_create_psi_mats(num_xi, PHI_mats,
    components, x_cat, C_vals, C_vecs, C_ranks, C_sizes, C_rsp);
}

parameters {
#include parameters.stan
  vector[sum(num_xi)] xi;
}

transformed parameters {
 vector[num_obs] f_latent[num_comps] = STAN_build_f(components,
  num_xi, C_ranks, seq_B, L, PSI_mats, alpha, ell, xi);
}

model {
  vector[num_obs] f_sum = STAN_vectorsum(f_latent, num_obs) + c_hat;
  xi ~ normal(0, 1);
#include model.stan
}
