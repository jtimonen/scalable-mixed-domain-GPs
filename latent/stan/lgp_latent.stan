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
#include functions-kernels.stan
}

data {
#include data.stan
  real<lower=0> delta;
}

transformed data{
  // Precompute fixed kernel matrices
  vector[num_obs] m0 = rep_vector(0.0, num_obs);
  matrix[num_obs, num_obs] K_const[num_comps] = STAN_kernel_const_all(
    num_obs, num_obs, x_cat, x_cat, x_cat_num_levels, components
  );
  vector[num_obs] delta_vec = rep_vector(delta, num_obs);
}

parameters {
#include parameters.stan
  vector[num_obs] eta[num_comps]; // isotropic versions of func components
}

transformed parameters {
  vector[num_obs] f_latent[num_comps];
  {
    matrix[num_obs, num_obs] Delta = diag_matrix(delta_vec);
    matrix[num_obs, num_obs] KX[num_comps] = STAN_kernel_all(
      num_obs, num_obs, K_const, components, x_cont, x_cont, alpha, ell
    );
    for(j in 1:num_comps){
      f_latent[j] = cholesky_decompose(KX[j] + Delta) * eta[j];
    }
  }
}

model {
  vector[num_obs] f_sum = STAN_vectorsum(f_latent, num_obs) + c_hat;
  for(j in 1:num_comps){ target += normal_lpdf(eta[j] | 0, 1);}
#include model.stan
}
