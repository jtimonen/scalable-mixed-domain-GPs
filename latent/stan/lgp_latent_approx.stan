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
  int idx_z[num_comps] = components[:, 8];
  matrix[num_obs, num_bf] PHI_mats[num_cov_cont] = 
    STAN_create_phi_mats(num_obs, num_bf, num_cov_cont, scale_bf, x_cont, X_hr);
  matrix[num_obs, sum(num_xi)] PSI = STAN_create_psi_mats();
  vector[num_bf*sum(C_ranks)] seq_C;
}

parameters {
#include parameters.stan
  vector[sum(num_xi)] xi;
}

transformed parameters {
  vector[num_obs] f_latent[num_comps];
  {
    int idx_ell = 0;

    // Build the components
    for(j in 1:num_comps) {
      vector[num_xi[j]] dj;
      if(components[j,1]==0) {
        d_j = rep_vector(alpha[j], C_ranks[j]);
      } else {
        idx_ell += 1;
        int ix = components[j, 9];
        real L = X_hr[ix] * scale_bf;
        dj = STAN_diag_spd_eq(alpha[j], ell[idx_ell], seq_Cj, Lj)
      }
      f_latent[j] = PSI_j * (dj .* STAN_subvector(xi, num_xi, j));
    }
  }
}

model {
  vector[num_obs] f_sum = STAN_vectorsum(f_latent, num_obs) + c_hat;
  xi ~ normal(0, 1);
#include model.stan
}
