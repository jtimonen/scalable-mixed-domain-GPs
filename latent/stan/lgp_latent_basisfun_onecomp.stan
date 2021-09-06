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
}

transformed data{
  matrix[num_obs, num_bf] PHI_mats[num_cov_cont];
  vector[num_bf] seq_C1 = STAN_seq_len_rep_times(num_bf, 1);
  real L = 1.0 * scale_bf;
  for(ix in 1:num_cov_cont) {
    PHI_mats[ix] = STAN_PHI_eq(x_cont[ix], seq_C1, L);
  }
}

parameters {
#include parameters.stan
  vector[num_bf] xi_1[1];
}

transformed parameters {
  vector[num_obs] f_latent[num_comps];
  {
    vector[num_bf] d1 = STAN_diag_spd_eq(alpha[1], ell[1], seq_C1, L);
    f_latent[1] = PHI_mats[1] * (d1 .* xi_1[1]);   // (N,M) x (M)
  }
}

model {
  vector[num_obs] f_sum = STAN_vectorsum(f_latent, num_obs) + c_hat;
  xi_1[1] ~ normal(0, 1);
#include model.stan
}
