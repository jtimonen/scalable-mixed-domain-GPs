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
#include data.stan
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
  matrix[num_obs, 1*C3_num] PSI_3;
  vector[num_bf] seq_C1 = STAN_seq_len_rep_times(num_bf, 1);
  vector[num_bf*C2_num] seq_C2 = STAN_seq_len_rep_times(num_bf, C2_num);
  //vector[0*C3_num] seq_C3 = STAN_seq_len_rep_times(0, C3_num);
  
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
    int i1 = 1 + (c-1)*1;
    int i2 = c*1;
    vector[num_obs] vpc = sqrt(C3_vals[c]) * C3_vecs[x_cat[idx_z3, :],c];
    PSI_3[:, i1:i2] = rep_matrix(vpc, 1);
  }
  
}

parameters {
#include parameters.stan
  vector[num_bf] xi_1;
  vector[num_bf*C2_num] xi_2;
  vector[C3_num] xi_3;
}

transformed parameters {
  vector[num_obs] f_latent[num_comps];
  {
    
    // Multipliers
    vector[num_bf] d1 = STAN_diag_spd_eq(alpha[1], ell[1], seq_C1, L);
    vector[num_bf*C2_num] d2 = STAN_diag_spd_eq(alpha[2], ell[2], seq_C2, L);
    vector[C3_num] d3 = rep_vector(alpha[3], C3_num); //STAN_diag_spd_eq(alpha[3], ell[3], seq_C3, L);

    // Build the three components
    f_latent[1] = PSI_1 * (d1 .* xi_1);
    f_latent[2] = PSI_2 * (d2 .* xi_2);
    f_latent[3] = PSI_3 * (d3 .* xi_3);
  }
}

model {
  vector[num_obs] f_sum = STAN_vectorsum(f_latent, num_obs) + c_hat;
  
  // Multiplier priors
  xi_1 ~ normal(0, 1);
  xi_2 ~ normal(0, 1);
  xi_3 ~ normal(0, 1);
  
#include model.stan
}
