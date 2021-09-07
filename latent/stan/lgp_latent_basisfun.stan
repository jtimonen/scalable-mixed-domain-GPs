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


  // Get jth subvector
  // - vecs contains all vectors concatenated
  // - lens contains lengths of the vectors
  vector STAN_subvec(vector vecs, int[] lens, int j) {
    int i1 = sum(lens[1:(j-1)]) + 1;
    int i2 = sum(lens[1:j]);
    return vecs[i1:i2];
  }
  
  // Get ith element of jth vector
  // - vecs contains all vectors concatenated
  // - lens contains lengths of the vectors
  real STAN_elem_i_of_vec_j(real[] vecs, int[] lens, int j, int i) {
    if ((i > lens[j]) || (i <= 0)) {
      reject("element index <", i, "> out of bounds for vector <", j, ">");
    }
    return vecs[sum(lens[1:(j-1)]) + i];
  }
  
  // Get ith column of jth matrix
  // - mats contains all matrices flattened along columns and concatenated
  // - num_elems contains numbers of elements in each matrix
  // - num_rows contains numbers of rows in each matrix
  vector STAN_col_i_of_mat_j(vector mats, int[] num_elems, int[] num_rows,
  int j, int i) {
    int num_cols = num_elems[j] / num_rows[j];
    int idx = sum(num_elems[1:(j-1)]) + num_rows[j]*(i-1);
    if(i > num_cols || i <= 0) {
      reject("column index <", i, "> out of bounds for matrix <", j, ">");
    }
    return mats[(idx+1):(idx+num_rows[j])];
  }
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
  matrix[num_obs, num_bf] PHI_mats[num_cov_cont];
  matrix[num_obs, num_bf*sum(C_ranks)] PSI;
  vector[num_bf*sum(C_ranks)] seq_C;
  for(j in 1:num_comps) {
    int i1 = sum(C_ranks[1:(j-1)])*NUM_BF[j] + 1;
    int i2 = sum(C_ranks[1:j])*NUM_BF[j];
    seq_C[]
  }
  vector[num_bf*C2_rank] seq_C2 = STAN_seq_len_rep_times(num_bf, C2_rank);
  //vector[0*C3_num] seq_C3 = STAN_seq_len_rep_times(0, C3_num);
  
  // PHI
  for(ix in 1:num_cov_cont) {
    real L = X_hr[ix] * scale_bf;
    PHI_mats[ix] = STAN_PHI_eq(x_cont[ix], seq_C1, L);
  }
  
  // PSI_1
  PSI_1 = PHI_mats[1];
  
  // PSI_2
  for(c in 1:C2_rank) {
    int i1 = 1 + (c-1)*num_bf;
    int i2 = c*num_bf;
    vector[num_obs] vpc = sqrt(C2_vals[c]) * C2_vecs[c][x_cat[idx_z2, :]];
    PSI_2[:, i1:i2] = PHI_mats[1] .* rep_matrix(vpc, num_bf);
  }
  
    // PSI_3
  for(c in 1:C3_rank) {
    int i1 = 1 + (c-1)*1;
    int i2 = c*1;
    vector[num_obs] vpc = sqrt(C3_vals[c]) * C3_vecs[c][x_cat[idx_z3, :]];
    PSI_3[:, i1:i2] = rep_matrix(vpc, 1);
  }
  
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
