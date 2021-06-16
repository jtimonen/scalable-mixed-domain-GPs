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

  // Log prior density to be added to target
  real STAN_log_prior(real x, data int[] types, data real[] p) {
    real log_prior = 0;
    real t = x;
  
    // Possible transform and log of its absolute derivative
    if (types[2]==1){
      log_prior += log(fabs(2*x));
      t = square(x);
    }
    
    // Value of pdf
    if (types[1]==2){
      log_prior += normal_lpdf(t | p[1], p[2]); // 2 = normal
    }else if (types[1]==3){
      log_prior += student_t_lpdf(t | p[1], 0.0, 1.0); // 3 = student-t
    }else if (types[1]==4){
      log_prior += gamma_lpdf(t | p[1], p[2]); // 4 = gamma
    }else if (types[1]==5){
      log_prior += inv_gamma_lpdf(t | p[1], p[2]); // 5 = inv-gamma
    }else if (types[1]==6){
      log_prior += lognormal_lpdf(t | p[1], p[2]); // 6 = log-normal
    }
    
    return(log_prior);
  }
  
  // Categorical zero-sum kernel
  matrix STAN_kernel_zerosum(data int[] x1, data int[] x2, data int num_cat) {
    int n1 = size(x1); 
    int n2 = size(x2);
    matrix[n1, n2] K;
    for (i in 1:n1) {
      for (j in 1:n2) {
        if (x1[i] == x2[j]) {
          K[i,j] = 1;
        } else {
          K[i,j] = - inv(num_cat - 1); 
        }
      }
    }
    return(K);
  }
  
  // Binary mask kernel
  matrix STAN_kernel_cat(data int[] x1, data int[] x2) {
    int n1 = size(x1);
    int n2 = size(x2);
    matrix[n1,n2] K;
    for (i in 1:n1) {
      for (j in 1:n2) {
        K[i,j] = (x1[i] == x2[j]);
      }
    }
    return(K);
  }
  
  // BINARY MASK KERNEL
  matrix STAN_kernel_bin(data int[] x1, data int[] x2) {
    int n1 = size(x1);
    int n2 = size(x2);
    matrix[n1,n2] K;
    for (i in 1:n1) {
      for (j in 1:n2) {
        K[i,j] = (x1[i] == 0) * (x2[j] == 0);
      }
    }
    return(K);
  }
  
  // Compute one constant kernel matrix. Does not depend on parameters and
  // therefore this function never needs to be evaluated during sampling.
  matrix STAN_kernel_const(data int[] x1, data int[] x2, 
    data int kernel_type,  data int ncat) 
  {
    int n1 = num_elements(x1);
    int n2 = num_elements(x2);
    matrix[n1, n2] K;
    if (kernel_type == 1) {
      K = STAN_kernel_cat(x1, x2);
    } else if (kernel_type == 2) {
      K = STAN_kernel_bin(x1, x2);
    } else {
      // kernel_type = 0
      K = STAN_kernel_zerosum(x1, x2, ncat);
    }
    return(K);
  }
  
  // Compute all constant kernel matrices. These do not depend on parameters and
  // therefore this function never needs to be evaluated during sampling.
  matrix[] STAN_kernel_const_all(
    data int n1,           
    data int n2,
    data int[,] x1,        data int[,] x2,
    data int[] num_levels, data int[,] components)
  {
    int num_comps = size(components);
    matrix[n1, n2] K_const[num_comps];
    for (j in 1:num_comps) {
      matrix[n1, n2] K;
      int opts[9] = components[j];
      int ctype = opts[1];
      int ktype = opts[2];
      int idx_cat = opts[8];
      int idx_cont = opts[9];
      K = rep_matrix(1.0, n1, n2);
      
      // Compute kernel for categorical covariate
      if (ctype == 0 || ctype == 2) {
        int M = num_levels[idx_cat];
        K = K .* STAN_kernel_const(x1[idx_cat], x2[idx_cat], ktype, M);
      }
      K_const[j] = K;
    }
    return(K_const);
  }
  
  // Exponentiated quadratic kernel (with vector inputs)
  matrix STAN_kernel_eq(vector x1, vector x2, real alpha, real ell) {
    return(cov_exp_quad(to_array_1d(x1), to_array_1d(x2), alpha, ell));
  }
  
  /* 
    Compute all kernel matrices. These depend on parameters and
    therefore this function needs to be evaluated repeatedly during sampling.
  */
  matrix[] STAN_kernel_all(
    data int n1,
    data int n2,
    data matrix[] K_const,
    data int[,] components,
    data vector[] x1,
    data vector[] x2,
    real[] alpha,
    real[] ell)
  {
    int idx_ell = 0;
    int idx_alpha = 0;
    int num_comps = size(components);
    matrix[n1, n2] KX[num_comps];
  
    // Loop through components
    for(j in 1:num_comps){
      
      // 1. Initialize with constant part of the kernel matrix
      matrix[n1, n2] K = K_const[j];
      vector[n1] X1;
      vector[n2] X2;
  
      // 2. Get component properties
      int opts[9] = components[j];
      int ctype = opts[1];
      int idx_cont = opts[9];
      
      // 3. Pick the possible continuous covariate of this component
      if(ctype != 0){
        X1 = x1[idx_cont];
        X2 = x2[idx_cont];
      }
      
      // Compute the kernel matrix
      idx_alpha += 1;
      if(ctype != 0){
        idx_ell += 1;
        K = K .* STAN_kernel_eq(X1, X2, alpha[idx_alpha], ell[idx_ell]);
      } else {
        K = square(alpha[idx_alpha]) * K;
      }
      KX[j] = K; // store kernel matrix
    }
    return(KX);
  }
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
  // Precompute fixed kernel matrices
  vector[num_obs] m0 = rep_vector(0.0, num_obs);
  matrix[num_obs, num_obs] K_const[num_comps] = STAN_kernel_const_all(
    num_obs, num_obs, x_cat, x_cat, x_cat_num_levels, components
  );
  vector[num_obs] delta_vec = rep_vector(delta, num_obs);
  
  // Compute L for each continuous covariate
  real L_bf[num_cov_cont];
  for(j in 1:num_cov_cont) {
    L_bf[j] = c_bf*max(x_cont[j]);
  }
}

parameters {
  real<lower=1e-12> alpha[num_comps]; // component magnitudes
  real<lower=1e-12> ell[num_ell]; // lengthscales
  real<lower=1e-12> sigma[1]; // noise std
}

model {

  // Likelihood
  matrix[num_obs, num_obs] Ky = diag_matrix(delta_vec);
  matrix[num_obs, num_obs] KX[num_comps] = STAN_kernel_all(num_obs, num_obs,
      K_const, components, x_cont, x_cont, alpha, ell);
  for(j in 1:num_comps){
    Ky += KX[j];
  }
  Ky = add_diag(Ky, square(sigma[1]));
  y_norm ~ multi_normal_cholesky(m0, cholesky_decompose(Ky));
  
  // Priors
  for(j in 1:num_comps){
    target += STAN_log_prior(alpha[j], prior_alpha[j], hyper_alpha[j]);
  }
  for(j in 1:num_ell){
    target += STAN_log_prior(ell[j], prior_ell[j], hyper_ell[j]);
  }
  target += STAN_log_prior(sigma[1], prior_sigma[1], hyper_sigma[1]);
}

