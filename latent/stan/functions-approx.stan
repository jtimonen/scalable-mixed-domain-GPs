  // Create vector with elements 1, ..., N
  vector STAN_seq_len(data int N) {
    vector[N] v = rep_vector(1.0, N);
    for(n in 2:N) v[n] = n;
    return(v);
  }
  
  // Same as rep(x, times=N) in R
  vector STAN_rep_vector_times(data vector x, data int N) {
    return to_vector(rep_matrix(x, N));
  }
  
  // Get jth subvector
  // - vecs contains all vectors concatenated
  // - lens contains lengths of the vectors
  vector STAN_subvec(vector vecs, data int[] lens, int j) {
    int i1 = sum(lens[1:(j-1)]) + 1;
    int i2 = sum(lens[1:j]);
    return vecs[i1:i2];
  }
  
  // Get jth submatrix
  // - mats contains all matrices concatenated along second dimension
  // - num_cols contains numbers of columns of the matrices
  matrix STAN_submat(matrix mats, data int[] num_cols, int j) {
    int i1 = sum(num_cols[1:(j-1)]) + 1;
    int i2 = sum(num_cols[1:j]);
    return mats[:, i1:i2];
  }
  
  // Get ith element of jth vector
  // - vecs contains all vectors concatenated
  // - lens contains lengths of the vectors
  real STAN_elem_i_of_vec_j(real[] vecs, data int[] lens, int j, int i) {
    if ((i > lens[j]) || (i <= 0)) {
      reject("element index <", i, "> out of bounds for vector <", j, ">");
    }
    return vecs[sum(lens[1:(j-1)]) + i];
  }
  
  // Get ith column of jth matrix
  // - mats contains all matrices flattened along columns and concatenated
  // - num_elems contains numbers of elements in each matrix
  // - num_rows contains numbers of rows in each matrix
  vector STAN_col_i_of_mat_j(vector mats, data int[] num_elems, 
      data int[] num_rows, int j, int i) 
  {
    real num_cols = (1.0 * num_elems[j]) / (1.0 * num_rows[j]);
    int idx = sum(num_elems[1:(j-1)]) + num_rows[j]*(i-1);
    if(i > num_cols || i <= 0) {
      reject("column index <", i, "> out of bounds for matrix <", j, ">");
    }
    return mats[(idx+1):(idx+num_rows[j])];
  }
  
  // Basis function matrix (EQ kernel)
  matrix STAN_basisfun_eq(vector x, data matrix mat_B, data real L) {
    int N = num_elements(x);
    int B = cols(mat_B);
    matrix[N, B] mat_X = rep_matrix(x+L, B);
    return 1.0/sqrt(L)*sin(0.5*pi()/L * mat_X .* mat_B);
  }
  
  // Compute spectral density of EQ kernel
  vector STAN_log_spd_eq(real alpha, real ell, vector omega){
    return 2*log(alpha)+log(ell)+0.5*log(2*pi())-0.5*ell^2*omega .* omega;
  }
  
  // Compute the multipliers s_b
  vector STAN_basisfun_eq_multipliers(real alpha, real ell, data vector seq_B,
      real L)
  {
    return sqrt(1e-20+exp(STAN_log_spd_eq(alpha, ell, 0.5*pi()*seq_B/L)));
  }
  
  // Create all PHI matrices
  matrix[] STAN_create_basisfun_mats(data vector[] x_cont, data matrix mat_B,
    data vector L) 
  {
    int N = rows(mat_B);
    int B = cols(mat_B);
    int D = size(x_cont);
    matrix[N, B] PHI_mats[D];
    for(ix in 1:D) {
      PHI_mats[ix] = STAN_basisfun_eq(x_cont[ix], mat_B, L[ix]);
    }
    return(PHI_mats);
  }
  
  // Create psi matrix with component type 0
  matrix STAN_create_psi_type0(int N, int j, data int[] z, 
      data real[] C_vals, data vector C_vecs, data int[] C_ranks,
      data int[] C_sizes, data int[] C_rsp) 
  {
    int R = C_ranks[j];
    int S = C_sizes[j];
    matrix[N, R] PSI_j;
    for(r in 1:R) {
      real eval = STAN_elem_i_of_vec_j(C_vals, C_ranks, j, r);
      vector[S] evec = STAN_col_i_of_mat_j(C_vecs, C_rsp, C_sizes, j, r);
      PSI_j[:,r] = sqrt(eval) * evec[z];
    }
    return(PSI_j);
  }
  
  // Create psi matrix with component type 2
  matrix STAN_create_psi_type2(int N, int j, data int[] z, 
      data real[] C_vals, data vector C_vecs, data int[] C_ranks,
      data int[] C_sizes, data int[] C_rsp, int M, matrix PHI) 
  {
    int R = C_ranks[j];
    int S = C_sizes[j];
    matrix[N, M*R] PSI_j;
    for(r in 1:R) {
      real eval = STAN_elem_i_of_vec_j(C_vals, C_ranks, j, r);
      vector[S] evec = STAN_col_i_of_mat_j(C_vecs, C_rsp, C_sizes, j, r);
      int j1 = 1 + (r-1)*M;
      int j2 = r*M;
      PSI_j[:,j1:j2] = PHI .* rep_matrix(sqrt(eval) * evec[z], M);
    }
    return(PSI_j);
  }
  
  // Create all PSI matrices
  // - returns a matrix of size N x total_num_xi
  matrix STAN_create_psi_mats(data int[] num_xi, data matrix[] PHI_mats,
      data int[,] components, data int[,] x_cat, data real[] C_vals,
      data vector C_vecs, data int[] C_ranks, data int[] C_sizes, 
      data int[] C_rsp) 
  {
    int N = rows(PHI_mats[1]);
    int B = cols(PHI_mats[1]);
    matrix[N, sum(num_xi)] PSI;
    int J = size(components);
    for(j in 1:J) {
      int ctype = components[j,1];
      int idx_x = components[j,9];
      int idx_z = components[j,8];
      int i1 = sum(num_xi[1:(j-1)]) + 1;
      int i2 = sum(num_xi[1:j]);
      matrix[N, num_xi[j]] PSI_j;
      if (ctype==0)  {
        PSI_j = STAN_create_psi_type0(N, j, x_cat[idx_z, :], C_vals, C_vecs,
          C_ranks, C_sizes, C_rsp);
      } else if (ctype==1) {
        PSI_j = PHI_mats[idx_x];
      } else {
        PSI_j = STAN_create_psi_type2(N, j, x_cat[idx_z, :], C_vals, C_vecs,
          C_ranks, C_sizes, C_rsp, B, PHI_mats[idx_x]);
      }
      PSI[:, i1:i2] = PSI_j;
    }
    return(PSI);
  }
  
  // Build the components of the latent signal f
  vector[] STAN_build_f_latent(data int[,] components, 
      data int[] num_xi, data int[] C_ranks, data vector seq_B, 
      data vector L, data matrix PSI, real[] alpha, real[] ell, vector xi)
  {
    int J = size(components);
    int N = rows(PSI);
    int B = num_elements(seq_B);
    vector[N] f_latent[J];
    int idx_ell = 0;

    // Build the components
    for(j in 1:J) {
      int ctype = components[j, 1];
      int idx_x = components[j, 9];
      int R = C_ranks[j];
      vector[num_xi[j]] dj;
      if(ctype==0) {
        dj = rep_vector(alpha[j], R);
      } else {
        idx_ell += 1;
        dj = STAN_basisfun_eq_multipliers(alpha[j], ell[idx_ell], STAN_rep_vector_times(seq_B, R), L[idx_x]);
      }
      f_latent[j] = STAN_submat(PSI, num_xi, j) * (dj .* STAN_subvec(xi, num_xi, j));
    }
    return(f_latent);
  }
