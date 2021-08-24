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
