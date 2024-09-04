functions {
  // --------------- utils ---------------
  
  // Create vector with elements 1, ..., N
  vector seq_len(data int N) {
    vector[N] v = rep_vector(1.0, N);
    for (n in 2 : N) 
      v[n] = n;
    return (v);
  }
  
  // Create integer array with elements 1, ..., N
  array[] int seq_len_int(data int N) {
    array[N] int v = rep_array(1, N);
    for (n in 2 : N) 
      v[n] = n;
    return (v);
  }
  
  // Same as rep(x, times=N) in R
  vector rep_vector_times(vector x, data int N) {
    return to_vector(rep_matrix(x, N));
  }
  
  // Repeat a matrix horizontally
  matrix rep_matrix_times(data matrix x, data int times) {
    int N = rows(x);
    int M = cols(x);
    matrix[N, M * times] out;
    for (j in 1 : times) {
      out[ : , (1 + (j - 1) * M) : (j * M)] = x;
    }
    return (out);
  }
  
  // Find indices of val
  array[] int which(array[] int vals, int val) {
    int n = size(vals);
    array[n] int out;
    int J = 0;
    for (i in 1 : n) {
      if (i == val) {
        out[n] = i;
        J = J + 1;
      }
    }
    return (out[1 : J]);
  }
  
  // Find first occurence of value in array
  int first(array[] int arr, int val) {
    int J = size(arr);
    for (j in 1 : J) {
      if (arr[j] == val) {
        return (j);
      }
    }
    return (0);
  }
  
  // Find last occurence of value in array
  int last(array[] int arr, int val) {
    int J = size(arr);
    array[J] int r_arr = reverse(arr);
    for (j in 1 : J) {
      if (r_arr[j] == val) {
        return (J - j + 1);
      }
    }
    return (0);
  }
  
  // Create an index array where first column is start indices and second has
  // the end indices
  // Assumes that arr is sorted
  array[,] int inds_array(array[] int arr, int G) {
    array[G, 2] int out;
    for (g in 1 : G) {
      out[g, 1] = first(arr, g);
      out[g, 2] = last(arr, g);
    }
    print("Created index array:");
    return (out);
  }
  
  // Transform parameter to log natural scale (vectorized)
  // - all vector arguments should have the same length
  vector to_log_natural_scale(vector log_z, vector log_mu, vector log_sigma) {
    return (log_sum_exp(log_z + log_sigma, log_mu));
  }
  
  // --------------- gp/eigenfunctions ---------------
  // Eigen functions for shared term
  matrix eigfun_shared(vector x, data matrix mat_B, data real L) {
    return (bf_eq(x, mat_B, L));
  }
  
  // Eigen functions for grouped term
  matrix eigfun_grouped(vector x, data matrix mat_B, data real L,
                        data array[] int z, data int G) {
    int N = num_elements(x);
    int B = cols(mat_B);
    matrix[N, B] PHI = bf_eq(x, mat_B, L);
    matrix[N, G - 1] VARPHI = bf_zs(z, G);
    return (bf_eqXzs(PHI, VARPHI));
  }
  
  // Basis function matrix (EQ kernel)
  matrix bf_eq(vector x, data matrix mat_B, data real L) {
    int N = num_elements(x);
    int B = cols(mat_B);
    matrix[N, B] mat_X = rep_matrix(x + L, B);
    matrix[N, B] PHI = 1.0 / sqrt(L) * sin(0.5 * pi() / L * mat_X .* mat_B);
    return (PHI);
  }
  
  // Basis function matrix (zerosum kernel)
  matrix bf_zs(data array[] int z, data int G) {
    int N = size(z);
    int Gm1 = G - 1;
    matrix[N, G - 1] VARPHI;
    matrix[G, G - 1] H = helmert_norm(G);
    for (g in 2 : G) {
      VARPHI[ : , g - 1] = H[z, g - 1];
    }
    return (VARPHI);
  }
  
  // Compute product basis functions (interaction)
  // - PHI = the evaluated basis functions (size num_obs x num_bf)
  // - VARPHI = the evaluated cat basis functions (size num_obs x (G-1))
  matrix bf_eqXzs(matrix PHI, matrix VARPHI) {
    int N = rows(PHI);
    int B = cols(PHI);
    int Gm1 = cols(VARPHI);
    matrix[N, B * Gm1] PSI;
    int idx;
    for (g in 1 : Gm1) {
      for (b in 1 : B) {
        idx = (g - 1) * B + b;
        PSI[ : , idx] = PHI[ : , b] .* VARPHI[ : , g];
      }
    }
    return (PSI);
  }
  
  // Helmert contrast matrix
  matrix helmert(int G) {
    matrix[G, G - 1] H;
    for (r in 1 : G) {
      for (c in 2 : G) {
        if (r < c) {
          H[r, c - 1] = -1;
        } else if (r == c) {
          H[r, c - 1] = c - 1;
        } else {
          H[r, c - 1] = 0;
        }
      }
    }
    return (H);
  }
  
  // Helmert contrast matrix with normalized columns
  matrix helmert_norm(int G) {
    matrix[G, G - 1] H = helmert(G);
    for (g in 2 : G) {
      H[ : , g - 1] = H[ : , g - 1] / norm2(H[ : , g - 1]);
    }
    return (H);
  }
  
  // --------------- gp/eigenvalues ---------------
  // Square roots of Eigenvalues
  vector eigval_shared(real alpha, real ell, data vector seq_B, data real L) {
    return (bf_eq_multips(alpha, ell, seq_B, L));
  }
  
  // Square roots of Eigenvalues
  vector eigval_grouped(real alpha, real ell, data vector seq_B, data real L,
                        data int G) {
    int B = num_elements(seq_B);
    int Gm1 = G - 1;
    vector[B] s = bf_eq_multips(alpha, ell, seq_B, L);
    vector[B * Gm1] DELTA;
    for (g in 1 : Gm1) {
      DELTA[((g - 1) * B + 1) : (g * B)] = sqrt(1.0 + 1.0 / Gm1) * s;
    }
    return (DELTA);
  }
  
  // Compute spectral density of EQ kernel
  vector log_spd_eq(real alpha, real ell, vector omega) {
    return 2 * log(alpha) + log(ell) + 0.5 * log(2 * pi())
           - 0.5 * ell ^ 2 * omega .* omega;
  }
  
  // Compute the multipliers s_b
  vector bf_eq_multips(real alpha, real ell, data vector seq_B, data real L) {
    return exp(0.5 * log_spd_eq(alpha, ell, 0.5 * pi() * seq_B / L));
  }
}
data {
  int<lower=0, upper=1> prior_only;
  
  // LonModel
  int<lower=1> n_LON; // number of input points (LON)
  
  // TermList for f_sum
  
  // Term f_baseline_id (GroupedOffsetTerm)
  array[n_LON] int dat_id_LON; // categ input
  int<lower=1> G_id; // num of groups
  
  // Term f_gp_age (GPTerm)
  vector[n_LON] dat_age_unit_LON; // continuous input
  real<lower=0> L_age; // domain size
  int<lower=1> B_age; // num of basis funs
  
  // Term f_gp_x (GPTerm)
  vector[n_LON] dat_x_unit_LON; // continuous input
  real<lower=0> L_x; // domain size
  int<lower=1> B_x; // num of basis funs
  
  // Term f_gp_x_u1 (GPTerm)
  vector[n_LON] dat_x_u1_unit_LON; // continuous input
  real<lower=0> L_x_u1; // domain size
  int<lower=1> B_x_u1; // num of basis funs
  
  // Term f_gp_x_u2 (GPTerm)
  vector[n_LON] dat_x_u2_unit_LON; // continuous input
  real<lower=0> L_x_u2; // domain size
  int<lower=1> B_x_u2; // num of basis funs
  
  // Term f_gp_x_u3 (GPTerm)
  vector[n_LON] dat_x_u3_unit_LON; // continuous input
  real<lower=0> L_x_u3; // domain size
  int<lower=1> B_x_u3; // num of basis funs
  
  // Term f_gp_x_u4 (GPTerm)
  vector[n_LON] dat_x_u4_unit_LON; // continuous input
  real<lower=0> L_x_u4; // domain size
  int<lower=1> B_x_u4; // num of basis funs
  
  // Term f_gp_x_u5 (GPTerm)
  vector[n_LON] dat_x_u5_unit_LON; // continuous input
  real<lower=0> L_x_u5; // domain size
  int<lower=1> B_x_u5; // num of basis funs
  
  // Term f_gp_x_u6 (GPTerm)
  vector[n_LON] dat_x_u6_unit_LON; // continuous input
  real<lower=0> L_x_u6; // domain size
  int<lower=1> B_x_u6; // num of basis funs
  
  // Term f_gp_x_u7 (GPTerm)
  vector[n_LON] dat_x_u7_unit_LON; // continuous input
  real<lower=0> L_x_u7; // domain size
  int<lower=1> B_x_u7; // num of basis funs
  
  // Term f_gp_x_u8 (GPTerm)
  vector[n_LON] dat_x_u8_unit_LON; // continuous input
  real<lower=0> L_x_u8; // domain size
  int<lower=1> B_x_u8; // num of basis funs
  
  // Term f_gp_x_u9 (GPTerm)
  vector[n_LON] dat_x_u9_unit_LON; // continuous input
  real<lower=0> L_x_u9; // domain size
  int<lower=1> B_x_u9; // num of basis funs
  
  // Term f_gp_x_u10 (GPTerm)
  vector[n_LON] dat_x_u10_unit_LON; // continuous input
  real<lower=0> L_x_u10; // domain size
  int<lower=1> B_x_u10; // num of basis funs
  
  // Term f_gp_x_u11 (GPTerm)
  vector[n_LON] dat_x_u11_unit_LON; // continuous input
  real<lower=0> L_x_u11; // domain size
  int<lower=1> B_x_u11; // num of basis funs
  
  // Term f_gp_x_u12 (GPTerm)
  vector[n_LON] dat_x_u12_unit_LON; // continuous input
  real<lower=0> L_x_u12; // domain size
  int<lower=1> B_x_u12; // num of basis funs
  
  // Term f_gp_x_u13 (GPTerm)
  vector[n_LON] dat_x_u13_unit_LON; // continuous input
  real<lower=0> L_x_u13; // domain size
  int<lower=1> B_x_u13; // num of basis funs
  
  // Term f_gp_x_u14 (GPTerm)
  vector[n_LON] dat_x_u14_unit_LON; // continuous input
  real<lower=0> L_x_u14; // domain size
  int<lower=1> B_x_u14; // num of basis funs
  
  // Term f_gp_x_u15 (GPTerm)
  vector[n_LON] dat_x_u15_unit_LON; // continuous input
  real<lower=0> L_x_u15; // domain size
  int<lower=1> B_x_u15; // num of basis funs
  
  // Term f_gp_x_u16 (GPTerm)
  vector[n_LON] dat_x_u16_unit_LON; // continuous input
  real<lower=0> L_x_u16; // domain size
  int<lower=1> B_x_u16; // num of basis funs
  
  // Term f_gp_ageXz (GPTerm)
  real<lower=0> L_ageXz; // domain size
  int<lower=1> B_ageXz; // num of basis funs
  array[n_LON] int dat_z_LON; // categ input
  int<lower=1> G_z; // num of groups
  
  // Term f_gp_ageXz_u1 (GPTerm)
  real<lower=0> L_ageXz_u1; // domain size
  int<lower=1> B_ageXz_u1; // num of basis funs
  array[n_LON] int dat_z_u1_LON; // categ input
  int<lower=1> G_z_u1; // num of groups
  
  // Term f_gp_ageXz_u2 (GPTerm)
  real<lower=0> L_ageXz_u2; // domain size
  int<lower=1> B_ageXz_u2; // num of basis funs
  array[n_LON] int dat_z_u2_LON; // categ input
  int<lower=1> G_z_u2; // num of groups
  
  // Term f_gp_ageXz_u3 (GPTerm)
  real<lower=0> L_ageXz_u3; // domain size
  int<lower=1> B_ageXz_u3; // num of basis funs
  array[n_LON] int dat_z_u3_LON; // categ input
  int<lower=1> G_z_u3; // num of groups
  
  // Term f_gp_ageXz_u4 (GPTerm)
  real<lower=0> L_ageXz_u4; // domain size
  int<lower=1> B_ageXz_u4; // num of basis funs
  array[n_LON] int dat_z_u4_LON; // categ input
  int<lower=1> G_z_u4; // num of groups
  
  // Term f_gp_ageXz_u5 (GPTerm)
  real<lower=0> L_ageXz_u5; // domain size
  int<lower=1> B_ageXz_u5; // num of basis funs
  array[n_LON] int dat_z_u5_LON; // categ input
  int<lower=1> G_z_u5; // num of groups
  
  // Term f_gp_ageXz_u6 (GPTerm)
  real<lower=0> L_ageXz_u6; // domain size
  int<lower=1> B_ageXz_u6; // num of basis funs
  array[n_LON] int dat_z_u6_LON; // categ input
  int<lower=1> G_z_u6; // num of groups
  
  // Term f_gp_ageXz_u7 (GPTerm)
  real<lower=0> L_ageXz_u7; // domain size
  int<lower=1> B_ageXz_u7; // num of basis funs
  array[n_LON] int dat_z_u7_LON; // categ input
  int<lower=1> G_z_u7; // num of groups
  
  // Term f_gp_ageXz_u8 (GPTerm)
  real<lower=0> L_ageXz_u8; // domain size
  int<lower=1> B_ageXz_u8; // num of basis funs
  array[n_LON] int dat_z_u8_LON; // categ input
  int<lower=1> G_z_u8; // num of groups
  
  // Term f_gp_ageXz_u9 (GPTerm)
  real<lower=0> L_ageXz_u9; // domain size
  int<lower=1> B_ageXz_u9; // num of basis funs
  array[n_LON] int dat_z_u9_LON; // categ input
  int<lower=1> G_z_u9; // num of groups
  
  // Term f_gp_ageXz_u10 (GPTerm)
  real<lower=0> L_ageXz_u10; // domain size
  int<lower=1> B_ageXz_u10; // num of basis funs
  array[n_LON] int dat_z_u10_LON; // categ input
  int<lower=1> G_z_u10; // num of groups
  
  // Term f_gp_ageXz_u11 (GPTerm)
  real<lower=0> L_ageXz_u11; // domain size
  int<lower=1> B_ageXz_u11; // num of basis funs
  array[n_LON] int dat_z_u11_LON; // categ input
  int<lower=1> G_z_u11; // num of groups
  
  // Term f_gp_ageXz_u12 (GPTerm)
  real<lower=0> L_ageXz_u12; // domain size
  int<lower=1> B_ageXz_u12; // num of basis funs
  array[n_LON] int dat_z_u12_LON; // categ input
  int<lower=1> G_z_u12; // num of groups
  
  // Term f_gp_ageXz_u13 (GPTerm)
  real<lower=0> L_ageXz_u13; // domain size
  int<lower=1> B_ageXz_u13; // num of basis funs
  array[n_LON] int dat_z_u13_LON; // categ input
  int<lower=1> G_z_u13; // num of groups
  
  // Term f_gp_ageXz_u14 (GPTerm)
  real<lower=0> L_ageXz_u14; // domain size
  int<lower=1> B_ageXz_u14; // num of basis funs
  array[n_LON] int dat_z_u14_LON; // categ input
  int<lower=1> G_z_u14; // num of groups
  
  // Term f_gp_ageXz_u15 (GPTerm)
  real<lower=0> L_ageXz_u15; // domain size
  int<lower=1> B_ageXz_u15; // num of basis funs
  array[n_LON] int dat_z_u15_LON; // categ input
  int<lower=1> G_z_u15; // num of groups
  
  // Term f_gp_ageXz_u16 (GPTerm)
  real<lower=0> L_ageXz_u16; // domain size
  int<lower=1> B_ageXz_u16; // num of basis funs
  array[n_LON] int dat_z_u16_LON; // categ input
  int<lower=1> G_z_u16; // num of groups
  
  // Observation model
  vector[n_LON] dat_y; // Longitudinal observations
  real c_hat; // fixed offset
}
transformed data {
  // LonModel
  
  // TermList for f_sum
  
  // Term f_gp_age (GPTerm)
  vector[B_age] seq_B_age = seq_len(B_age);
  matrix[n_LON, B_age] mat_B_age_LON = transpose(rep_matrix(seq_B_age, n_LON));
  matrix[n_LON, B_age] PSI_age_LON = eigfun_shared(dat_age_unit_LON,
                                                   mat_B_age_LON, L_age);
  
  // Term f_gp_x (GPTerm)
  vector[B_x] seq_B_x = seq_len(B_x);
  matrix[n_LON, B_x] mat_B_x_LON = transpose(rep_matrix(seq_B_x, n_LON));
  matrix[n_LON, B_x] PSI_x_LON = eigfun_shared(dat_x_unit_LON, mat_B_x_LON,
                                               L_x);
  
  // Term f_gp_x_u1 (GPTerm)
  vector[B_x_u1] seq_B_x_u1 = seq_len(B_x_u1);
  matrix[n_LON, B_x_u1] mat_B_x_u1_LON = transpose(rep_matrix(seq_B_x_u1,
                                                              n_LON));
  matrix[n_LON, B_x_u1] PSI_x_u1_LON = eigfun_shared(dat_x_u1_unit_LON,
                                                     mat_B_x_u1_LON, L_x_u1);
  
  // Term f_gp_x_u2 (GPTerm)
  vector[B_x_u2] seq_B_x_u2 = seq_len(B_x_u2);
  matrix[n_LON, B_x_u2] mat_B_x_u2_LON = transpose(rep_matrix(seq_B_x_u2,
                                                              n_LON));
  matrix[n_LON, B_x_u2] PSI_x_u2_LON = eigfun_shared(dat_x_u2_unit_LON,
                                                     mat_B_x_u2_LON, L_x_u2);
  
  // Term f_gp_x_u3 (GPTerm)
  vector[B_x_u3] seq_B_x_u3 = seq_len(B_x_u3);
  matrix[n_LON, B_x_u3] mat_B_x_u3_LON = transpose(rep_matrix(seq_B_x_u3,
                                                              n_LON));
  matrix[n_LON, B_x_u3] PSI_x_u3_LON = eigfun_shared(dat_x_u3_unit_LON,
                                                     mat_B_x_u3_LON, L_x_u3);
  
  // Term f_gp_x_u4 (GPTerm)
  vector[B_x_u4] seq_B_x_u4 = seq_len(B_x_u4);
  matrix[n_LON, B_x_u4] mat_B_x_u4_LON = transpose(rep_matrix(seq_B_x_u4,
                                                              n_LON));
  matrix[n_LON, B_x_u4] PSI_x_u4_LON = eigfun_shared(dat_x_u4_unit_LON,
                                                     mat_B_x_u4_LON, L_x_u4);
  
  // Term f_gp_x_u5 (GPTerm)
  vector[B_x_u5] seq_B_x_u5 = seq_len(B_x_u5);
  matrix[n_LON, B_x_u5] mat_B_x_u5_LON = transpose(rep_matrix(seq_B_x_u5,
                                                              n_LON));
  matrix[n_LON, B_x_u5] PSI_x_u5_LON = eigfun_shared(dat_x_u5_unit_LON,
                                                     mat_B_x_u5_LON, L_x_u5);
  
  // Term f_gp_x_u6 (GPTerm)
  vector[B_x_u6] seq_B_x_u6 = seq_len(B_x_u6);
  matrix[n_LON, B_x_u6] mat_B_x_u6_LON = transpose(rep_matrix(seq_B_x_u6,
                                                              n_LON));
  matrix[n_LON, B_x_u6] PSI_x_u6_LON = eigfun_shared(dat_x_u6_unit_LON,
                                                     mat_B_x_u6_LON, L_x_u6);
  
  // Term f_gp_x_u7 (GPTerm)
  vector[B_x_u7] seq_B_x_u7 = seq_len(B_x_u7);
  matrix[n_LON, B_x_u7] mat_B_x_u7_LON = transpose(rep_matrix(seq_B_x_u7,
                                                              n_LON));
  matrix[n_LON, B_x_u7] PSI_x_u7_LON = eigfun_shared(dat_x_u7_unit_LON,
                                                     mat_B_x_u7_LON, L_x_u7);
  
  // Term f_gp_x_u8 (GPTerm)
  vector[B_x_u8] seq_B_x_u8 = seq_len(B_x_u8);
  matrix[n_LON, B_x_u8] mat_B_x_u8_LON = transpose(rep_matrix(seq_B_x_u8,
                                                              n_LON));
  matrix[n_LON, B_x_u8] PSI_x_u8_LON = eigfun_shared(dat_x_u8_unit_LON,
                                                     mat_B_x_u8_LON, L_x_u8);
  
  // Term f_gp_x_u9 (GPTerm)
  vector[B_x_u9] seq_B_x_u9 = seq_len(B_x_u9);
  matrix[n_LON, B_x_u9] mat_B_x_u9_LON = transpose(rep_matrix(seq_B_x_u9,
                                                              n_LON));
  matrix[n_LON, B_x_u9] PSI_x_u9_LON = eigfun_shared(dat_x_u9_unit_LON,
                                                     mat_B_x_u9_LON, L_x_u9);
  
  // Term f_gp_x_u10 (GPTerm)
  vector[B_x_u10] seq_B_x_u10 = seq_len(B_x_u10);
  matrix[n_LON, B_x_u10] mat_B_x_u10_LON = transpose(rep_matrix(seq_B_x_u10,
                                                                n_LON));
  matrix[n_LON, B_x_u10] PSI_x_u10_LON = eigfun_shared(dat_x_u10_unit_LON,
                                                       mat_B_x_u10_LON,
                                                       L_x_u10);
  
  // Term f_gp_x_u11 (GPTerm)
  vector[B_x_u11] seq_B_x_u11 = seq_len(B_x_u11);
  matrix[n_LON, B_x_u11] mat_B_x_u11_LON = transpose(rep_matrix(seq_B_x_u11,
                                                                n_LON));
  matrix[n_LON, B_x_u11] PSI_x_u11_LON = eigfun_shared(dat_x_u11_unit_LON,
                                                       mat_B_x_u11_LON,
                                                       L_x_u11);
  
  // Term f_gp_x_u12 (GPTerm)
  vector[B_x_u12] seq_B_x_u12 = seq_len(B_x_u12);
  matrix[n_LON, B_x_u12] mat_B_x_u12_LON = transpose(rep_matrix(seq_B_x_u12,
                                                                n_LON));
  matrix[n_LON, B_x_u12] PSI_x_u12_LON = eigfun_shared(dat_x_u12_unit_LON,
                                                       mat_B_x_u12_LON,
                                                       L_x_u12);
  
  // Term f_gp_x_u13 (GPTerm)
  vector[B_x_u13] seq_B_x_u13 = seq_len(B_x_u13);
  matrix[n_LON, B_x_u13] mat_B_x_u13_LON = transpose(rep_matrix(seq_B_x_u13,
                                                                n_LON));
  matrix[n_LON, B_x_u13] PSI_x_u13_LON = eigfun_shared(dat_x_u13_unit_LON,
                                                       mat_B_x_u13_LON,
                                                       L_x_u13);
  
  // Term f_gp_x_u14 (GPTerm)
  vector[B_x_u14] seq_B_x_u14 = seq_len(B_x_u14);
  matrix[n_LON, B_x_u14] mat_B_x_u14_LON = transpose(rep_matrix(seq_B_x_u14,
                                                                n_LON));
  matrix[n_LON, B_x_u14] PSI_x_u14_LON = eigfun_shared(dat_x_u14_unit_LON,
                                                       mat_B_x_u14_LON,
                                                       L_x_u14);
  
  // Term f_gp_x_u15 (GPTerm)
  vector[B_x_u15] seq_B_x_u15 = seq_len(B_x_u15);
  matrix[n_LON, B_x_u15] mat_B_x_u15_LON = transpose(rep_matrix(seq_B_x_u15,
                                                                n_LON));
  matrix[n_LON, B_x_u15] PSI_x_u15_LON = eigfun_shared(dat_x_u15_unit_LON,
                                                       mat_B_x_u15_LON,
                                                       L_x_u15);
  
  // Term f_gp_x_u16 (GPTerm)
  vector[B_x_u16] seq_B_x_u16 = seq_len(B_x_u16);
  matrix[n_LON, B_x_u16] mat_B_x_u16_LON = transpose(rep_matrix(seq_B_x_u16,
                                                                n_LON));
  matrix[n_LON, B_x_u16] PSI_x_u16_LON = eigfun_shared(dat_x_u16_unit_LON,
                                                       mat_B_x_u16_LON,
                                                       L_x_u16);
  
  // Term f_gp_ageXz (GPTerm)
  vector[B_ageXz] seq_B_ageXz = seq_len(B_ageXz);
  matrix[n_LON, B_ageXz] mat_B_ageXz_LON = transpose(rep_matrix(seq_B_ageXz,
                                                                n_LON));
  matrix[n_LON, B_ageXz * (G_z - 1)] PSI_ageXz_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_LON,
                                                                    L_ageXz,
                                                                    dat_z_LON,
                                                                    G_z);
  
  // Term f_gp_ageXz_u1 (GPTerm)
  vector[B_ageXz_u1] seq_B_ageXz_u1 = seq_len(B_ageXz_u1);
  matrix[n_LON, B_ageXz_u1] mat_B_ageXz_u1_LON = transpose(rep_matrix(
                                                           seq_B_ageXz_u1,
                                                           n_LON));
  matrix[n_LON, B_ageXz_u1 * (G_z_u1 - 1)] PSI_ageXz_u1_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u1_LON,
                                                                    L_ageXz_u1,
                                                                    dat_z_u1_LON,
                                                                    G_z_u1);
  
  // Term f_gp_ageXz_u2 (GPTerm)
  vector[B_ageXz_u2] seq_B_ageXz_u2 = seq_len(B_ageXz_u2);
  matrix[n_LON, B_ageXz_u2] mat_B_ageXz_u2_LON = transpose(rep_matrix(
                                                           seq_B_ageXz_u2,
                                                           n_LON));
  matrix[n_LON, B_ageXz_u2 * (G_z_u2 - 1)] PSI_ageXz_u2_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u2_LON,
                                                                    L_ageXz_u2,
                                                                    dat_z_u2_LON,
                                                                    G_z_u2);
  
  // Term f_gp_ageXz_u3 (GPTerm)
  vector[B_ageXz_u3] seq_B_ageXz_u3 = seq_len(B_ageXz_u3);
  matrix[n_LON, B_ageXz_u3] mat_B_ageXz_u3_LON = transpose(rep_matrix(
                                                           seq_B_ageXz_u3,
                                                           n_LON));
  matrix[n_LON, B_ageXz_u3 * (G_z_u3 - 1)] PSI_ageXz_u3_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u3_LON,
                                                                    L_ageXz_u3,
                                                                    dat_z_u3_LON,
                                                                    G_z_u3);
  
  // Term f_gp_ageXz_u4 (GPTerm)
  vector[B_ageXz_u4] seq_B_ageXz_u4 = seq_len(B_ageXz_u4);
  matrix[n_LON, B_ageXz_u4] mat_B_ageXz_u4_LON = transpose(rep_matrix(
                                                           seq_B_ageXz_u4,
                                                           n_LON));
  matrix[n_LON, B_ageXz_u4 * (G_z_u4 - 1)] PSI_ageXz_u4_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u4_LON,
                                                                    L_ageXz_u4,
                                                                    dat_z_u4_LON,
                                                                    G_z_u4);
  
  // Term f_gp_ageXz_u5 (GPTerm)
  vector[B_ageXz_u5] seq_B_ageXz_u5 = seq_len(B_ageXz_u5);
  matrix[n_LON, B_ageXz_u5] mat_B_ageXz_u5_LON = transpose(rep_matrix(
                                                           seq_B_ageXz_u5,
                                                           n_LON));
  matrix[n_LON, B_ageXz_u5 * (G_z_u5 - 1)] PSI_ageXz_u5_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u5_LON,
                                                                    L_ageXz_u5,
                                                                    dat_z_u5_LON,
                                                                    G_z_u5);
  
  // Term f_gp_ageXz_u6 (GPTerm)
  vector[B_ageXz_u6] seq_B_ageXz_u6 = seq_len(B_ageXz_u6);
  matrix[n_LON, B_ageXz_u6] mat_B_ageXz_u6_LON = transpose(rep_matrix(
                                                           seq_B_ageXz_u6,
                                                           n_LON));
  matrix[n_LON, B_ageXz_u6 * (G_z_u6 - 1)] PSI_ageXz_u6_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u6_LON,
                                                                    L_ageXz_u6,
                                                                    dat_z_u6_LON,
                                                                    G_z_u6);
  
  // Term f_gp_ageXz_u7 (GPTerm)
  vector[B_ageXz_u7] seq_B_ageXz_u7 = seq_len(B_ageXz_u7);
  matrix[n_LON, B_ageXz_u7] mat_B_ageXz_u7_LON = transpose(rep_matrix(
                                                           seq_B_ageXz_u7,
                                                           n_LON));
  matrix[n_LON, B_ageXz_u7 * (G_z_u7 - 1)] PSI_ageXz_u7_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u7_LON,
                                                                    L_ageXz_u7,
                                                                    dat_z_u7_LON,
                                                                    G_z_u7);
  
  // Term f_gp_ageXz_u8 (GPTerm)
  vector[B_ageXz_u8] seq_B_ageXz_u8 = seq_len(B_ageXz_u8);
  matrix[n_LON, B_ageXz_u8] mat_B_ageXz_u8_LON = transpose(rep_matrix(
                                                           seq_B_ageXz_u8,
                                                           n_LON));
  matrix[n_LON, B_ageXz_u8 * (G_z_u8 - 1)] PSI_ageXz_u8_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u8_LON,
                                                                    L_ageXz_u8,
                                                                    dat_z_u8_LON,
                                                                    G_z_u8);
  
  // Term f_gp_ageXz_u9 (GPTerm)
  vector[B_ageXz_u9] seq_B_ageXz_u9 = seq_len(B_ageXz_u9);
  matrix[n_LON, B_ageXz_u9] mat_B_ageXz_u9_LON = transpose(rep_matrix(
                                                           seq_B_ageXz_u9,
                                                           n_LON));
  matrix[n_LON, B_ageXz_u9 * (G_z_u9 - 1)] PSI_ageXz_u9_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u9_LON,
                                                                    L_ageXz_u9,
                                                                    dat_z_u9_LON,
                                                                    G_z_u9);
  
  // Term f_gp_ageXz_u10 (GPTerm)
  vector[B_ageXz_u10] seq_B_ageXz_u10 = seq_len(B_ageXz_u10);
  matrix[n_LON, B_ageXz_u10] mat_B_ageXz_u10_LON = transpose(rep_matrix(
                                                             seq_B_ageXz_u10,
                                                             n_LON));
  matrix[n_LON, B_ageXz_u10 * (G_z_u10 - 1)] PSI_ageXz_u10_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u10_LON,
                                                                    L_ageXz_u10,
                                                                    dat_z_u10_LON,
                                                                    G_z_u10);
  
  // Term f_gp_ageXz_u11 (GPTerm)
  vector[B_ageXz_u11] seq_B_ageXz_u11 = seq_len(B_ageXz_u11);
  matrix[n_LON, B_ageXz_u11] mat_B_ageXz_u11_LON = transpose(rep_matrix(
                                                             seq_B_ageXz_u11,
                                                             n_LON));
  matrix[n_LON, B_ageXz_u11 * (G_z_u11 - 1)] PSI_ageXz_u11_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u11_LON,
                                                                    L_ageXz_u11,
                                                                    dat_z_u11_LON,
                                                                    G_z_u11);
  
  // Term f_gp_ageXz_u12 (GPTerm)
  vector[B_ageXz_u12] seq_B_ageXz_u12 = seq_len(B_ageXz_u12);
  matrix[n_LON, B_ageXz_u12] mat_B_ageXz_u12_LON = transpose(rep_matrix(
                                                             seq_B_ageXz_u12,
                                                             n_LON));
  matrix[n_LON, B_ageXz_u12 * (G_z_u12 - 1)] PSI_ageXz_u12_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u12_LON,
                                                                    L_ageXz_u12,
                                                                    dat_z_u12_LON,
                                                                    G_z_u12);
  
  // Term f_gp_ageXz_u13 (GPTerm)
  vector[B_ageXz_u13] seq_B_ageXz_u13 = seq_len(B_ageXz_u13);
  matrix[n_LON, B_ageXz_u13] mat_B_ageXz_u13_LON = transpose(rep_matrix(
                                                             seq_B_ageXz_u13,
                                                             n_LON));
  matrix[n_LON, B_ageXz_u13 * (G_z_u13 - 1)] PSI_ageXz_u13_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u13_LON,
                                                                    L_ageXz_u13,
                                                                    dat_z_u13_LON,
                                                                    G_z_u13);
  
  // Term f_gp_ageXz_u14 (GPTerm)
  vector[B_ageXz_u14] seq_B_ageXz_u14 = seq_len(B_ageXz_u14);
  matrix[n_LON, B_ageXz_u14] mat_B_ageXz_u14_LON = transpose(rep_matrix(
                                                             seq_B_ageXz_u14,
                                                             n_LON));
  matrix[n_LON, B_ageXz_u14 * (G_z_u14 - 1)] PSI_ageXz_u14_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u14_LON,
                                                                    L_ageXz_u14,
                                                                    dat_z_u14_LON,
                                                                    G_z_u14);
  
  // Term f_gp_ageXz_u15 (GPTerm)
  vector[B_ageXz_u15] seq_B_ageXz_u15 = seq_len(B_ageXz_u15);
  matrix[n_LON, B_ageXz_u15] mat_B_ageXz_u15_LON = transpose(rep_matrix(
                                                             seq_B_ageXz_u15,
                                                             n_LON));
  matrix[n_LON, B_ageXz_u15 * (G_z_u15 - 1)] PSI_ageXz_u15_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u15_LON,
                                                                    L_ageXz_u15,
                                                                    dat_z_u15_LON,
                                                                    G_z_u15);
  
  // Term f_gp_ageXz_u16 (GPTerm)
  vector[B_ageXz_u16] seq_B_ageXz_u16 = seq_len(B_ageXz_u16);
  matrix[n_LON, B_ageXz_u16] mat_B_ageXz_u16_LON = transpose(rep_matrix(
                                                             seq_B_ageXz_u16,
                                                             n_LON));
  matrix[n_LON, B_ageXz_u16 * (G_z_u16 - 1)] PSI_ageXz_u16_LON = eigfun_grouped(dat_age_unit_LON,
                                                                    mat_B_ageXz_u16_LON,
                                                                    L_ageXz_u16,
                                                                    dat_z_u16_LON,
                                                                    G_z_u16);
}
parameters {
  // LonModel
  real<lower=0> sigma; // noise std parameter 
  
  // TermList for f_sum
  // Term f_baseline_id (GroupedOffsetTerm)
  vector[G_id] offset_id;
  
  // Term f_gp_age (GPTerm)
  real<lower=1e-12> alpha_age; // magnitude
  real<lower=1e-12> ell_age; // lengthscale
  vector[B_age] xi_age; // auxiliary
  
  // Term f_gp_x (GPTerm)
  real<lower=1e-12> alpha_x; // magnitude
  real<lower=1e-12> ell_x; // lengthscale
  vector[B_x] xi_x; // auxiliary
  
  // Term f_gp_x_u1 (GPTerm)
  real<lower=1e-12> alpha_x_u1; // magnitude
  real<lower=1e-12> ell_x_u1; // lengthscale
  vector[B_x_u1] xi_x_u1; // auxiliary
  
  // Term f_gp_x_u2 (GPTerm)
  real<lower=1e-12> alpha_x_u2; // magnitude
  real<lower=1e-12> ell_x_u2; // lengthscale
  vector[B_x_u2] xi_x_u2; // auxiliary
  
  // Term f_gp_x_u3 (GPTerm)
  real<lower=1e-12> alpha_x_u3; // magnitude
  real<lower=1e-12> ell_x_u3; // lengthscale
  vector[B_x_u3] xi_x_u3; // auxiliary
  
  // Term f_gp_x_u4 (GPTerm)
  real<lower=1e-12> alpha_x_u4; // magnitude
  real<lower=1e-12> ell_x_u4; // lengthscale
  vector[B_x_u4] xi_x_u4; // auxiliary
  
  // Term f_gp_x_u5 (GPTerm)
  real<lower=1e-12> alpha_x_u5; // magnitude
  real<lower=1e-12> ell_x_u5; // lengthscale
  vector[B_x_u5] xi_x_u5; // auxiliary
  
  // Term f_gp_x_u6 (GPTerm)
  real<lower=1e-12> alpha_x_u6; // magnitude
  real<lower=1e-12> ell_x_u6; // lengthscale
  vector[B_x_u6] xi_x_u6; // auxiliary
  
  // Term f_gp_x_u7 (GPTerm)
  real<lower=1e-12> alpha_x_u7; // magnitude
  real<lower=1e-12> ell_x_u7; // lengthscale
  vector[B_x_u7] xi_x_u7; // auxiliary
  
  // Term f_gp_x_u8 (GPTerm)
  real<lower=1e-12> alpha_x_u8; // magnitude
  real<lower=1e-12> ell_x_u8; // lengthscale
  vector[B_x_u8] xi_x_u8; // auxiliary
  
  // Term f_gp_x_u9 (GPTerm)
  real<lower=1e-12> alpha_x_u9; // magnitude
  real<lower=1e-12> ell_x_u9; // lengthscale
  vector[B_x_u9] xi_x_u9; // auxiliary
  
  // Term f_gp_x_u10 (GPTerm)
  real<lower=1e-12> alpha_x_u10; // magnitude
  real<lower=1e-12> ell_x_u10; // lengthscale
  vector[B_x_u10] xi_x_u10; // auxiliary
  
  // Term f_gp_x_u11 (GPTerm)
  real<lower=1e-12> alpha_x_u11; // magnitude
  real<lower=1e-12> ell_x_u11; // lengthscale
  vector[B_x_u11] xi_x_u11; // auxiliary
  
  // Term f_gp_x_u12 (GPTerm)
  real<lower=1e-12> alpha_x_u12; // magnitude
  real<lower=1e-12> ell_x_u12; // lengthscale
  vector[B_x_u12] xi_x_u12; // auxiliary
  
  // Term f_gp_x_u13 (GPTerm)
  real<lower=1e-12> alpha_x_u13; // magnitude
  real<lower=1e-12> ell_x_u13; // lengthscale
  vector[B_x_u13] xi_x_u13; // auxiliary
  
  // Term f_gp_x_u14 (GPTerm)
  real<lower=1e-12> alpha_x_u14; // magnitude
  real<lower=1e-12> ell_x_u14; // lengthscale
  vector[B_x_u14] xi_x_u14; // auxiliary
  
  // Term f_gp_x_u15 (GPTerm)
  real<lower=1e-12> alpha_x_u15; // magnitude
  real<lower=1e-12> ell_x_u15; // lengthscale
  vector[B_x_u15] xi_x_u15; // auxiliary
  
  // Term f_gp_x_u16 (GPTerm)
  real<lower=1e-12> alpha_x_u16; // magnitude
  real<lower=1e-12> ell_x_u16; // lengthscale
  vector[B_x_u16] xi_x_u16; // auxiliary
  
  // Term f_gp_ageXz (GPTerm)
  real<lower=1e-12> alpha_ageXz; // magnitude
  real<lower=1e-12> ell_ageXz; // lengthscale
  vector[(G_z - 1) * B_ageXz] xi_ageXz;
  
  // Term f_gp_ageXz_u1 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u1; // magnitude
  real<lower=1e-12> ell_ageXz_u1; // lengthscale
  vector[(G_z_u1 - 1) * B_ageXz_u1] xi_ageXz_u1;
  
  // Term f_gp_ageXz_u2 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u2; // magnitude
  real<lower=1e-12> ell_ageXz_u2; // lengthscale
  vector[(G_z_u2 - 1) * B_ageXz_u2] xi_ageXz_u2;
  
  // Term f_gp_ageXz_u3 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u3; // magnitude
  real<lower=1e-12> ell_ageXz_u3; // lengthscale
  vector[(G_z_u3 - 1) * B_ageXz_u3] xi_ageXz_u3;
  
  // Term f_gp_ageXz_u4 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u4; // magnitude
  real<lower=1e-12> ell_ageXz_u4; // lengthscale
  vector[(G_z_u4 - 1) * B_ageXz_u4] xi_ageXz_u4;
  
  // Term f_gp_ageXz_u5 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u5; // magnitude
  real<lower=1e-12> ell_ageXz_u5; // lengthscale
  vector[(G_z_u5 - 1) * B_ageXz_u5] xi_ageXz_u5;
  
  // Term f_gp_ageXz_u6 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u6; // magnitude
  real<lower=1e-12> ell_ageXz_u6; // lengthscale
  vector[(G_z_u6 - 1) * B_ageXz_u6] xi_ageXz_u6;
  
  // Term f_gp_ageXz_u7 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u7; // magnitude
  real<lower=1e-12> ell_ageXz_u7; // lengthscale
  vector[(G_z_u7 - 1) * B_ageXz_u7] xi_ageXz_u7;
  
  // Term f_gp_ageXz_u8 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u8; // magnitude
  real<lower=1e-12> ell_ageXz_u8; // lengthscale
  vector[(G_z_u8 - 1) * B_ageXz_u8] xi_ageXz_u8;
  
  // Term f_gp_ageXz_u9 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u9; // magnitude
  real<lower=1e-12> ell_ageXz_u9; // lengthscale
  vector[(G_z_u9 - 1) * B_ageXz_u9] xi_ageXz_u9;
  
  // Term f_gp_ageXz_u10 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u10; // magnitude
  real<lower=1e-12> ell_ageXz_u10; // lengthscale
  vector[(G_z_u10 - 1) * B_ageXz_u10] xi_ageXz_u10;
  
  // Term f_gp_ageXz_u11 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u11; // magnitude
  real<lower=1e-12> ell_ageXz_u11; // lengthscale
  vector[(G_z_u11 - 1) * B_ageXz_u11] xi_ageXz_u11;
  
  // Term f_gp_ageXz_u12 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u12; // magnitude
  real<lower=1e-12> ell_ageXz_u12; // lengthscale
  vector[(G_z_u12 - 1) * B_ageXz_u12] xi_ageXz_u12;
  
  // Term f_gp_ageXz_u13 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u13; // magnitude
  real<lower=1e-12> ell_ageXz_u13; // lengthscale
  vector[(G_z_u13 - 1) * B_ageXz_u13] xi_ageXz_u13;
  
  // Term f_gp_ageXz_u14 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u14; // magnitude
  real<lower=1e-12> ell_ageXz_u14; // lengthscale
  vector[(G_z_u14 - 1) * B_ageXz_u14] xi_ageXz_u14;
  
  // Term f_gp_ageXz_u15 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u15; // magnitude
  real<lower=1e-12> ell_ageXz_u15; // lengthscale
  vector[(G_z_u15 - 1) * B_ageXz_u15] xi_ageXz_u15;
  
  // Term f_gp_ageXz_u16 (GPTerm)
  real<lower=1e-12> alpha_ageXz_u16; // magnitude
  real<lower=1e-12> ell_ageXz_u16; // lengthscale
  vector[(G_z_u16 - 1) * B_ageXz_u16] xi_ageXz_u16;
}
transformed parameters {
  // LonModel
  
  // TermList for f_sum
  // Term f_baseline_id (GroupedOffsetTerm)
  vector[n_LON] f_baseline_id_LON = offset_id[dat_id_LON];
  
  // Term f_gp_age (GPTerm)
  vector[B_age] DELTA_age = eigval_shared(alpha_age, ell_age, seq_B_age,
                                          L_age);
  vector[n_LON] f_gp_age_LON = PSI_age_LON * (xi_age .* DELTA_age);
  
  // Term f_gp_x (GPTerm)
  vector[B_x] DELTA_x = eigval_shared(alpha_x, ell_x, seq_B_x, L_x);
  vector[n_LON] f_gp_x_LON = PSI_x_LON * (xi_x .* DELTA_x);
  
  // Term f_gp_x_u1 (GPTerm)
  vector[B_x_u1] DELTA_x_u1 = eigval_shared(alpha_x_u1, ell_x_u1, seq_B_x_u1,
                                            L_x_u1);
  vector[n_LON] f_gp_x_u1_LON = PSI_x_u1_LON * (xi_x_u1 .* DELTA_x_u1);
  
  // Term f_gp_x_u2 (GPTerm)
  vector[B_x_u2] DELTA_x_u2 = eigval_shared(alpha_x_u2, ell_x_u2, seq_B_x_u2,
                                            L_x_u2);
  vector[n_LON] f_gp_x_u2_LON = PSI_x_u2_LON * (xi_x_u2 .* DELTA_x_u2);
  
  // Term f_gp_x_u3 (GPTerm)
  vector[B_x_u3] DELTA_x_u3 = eigval_shared(alpha_x_u3, ell_x_u3, seq_B_x_u3,
                                            L_x_u3);
  vector[n_LON] f_gp_x_u3_LON = PSI_x_u3_LON * (xi_x_u3 .* DELTA_x_u3);
  
  // Term f_gp_x_u4 (GPTerm)
  vector[B_x_u4] DELTA_x_u4 = eigval_shared(alpha_x_u4, ell_x_u4, seq_B_x_u4,
                                            L_x_u4);
  vector[n_LON] f_gp_x_u4_LON = PSI_x_u4_LON * (xi_x_u4 .* DELTA_x_u4);
  
  // Term f_gp_x_u5 (GPTerm)
  vector[B_x_u5] DELTA_x_u5 = eigval_shared(alpha_x_u5, ell_x_u5, seq_B_x_u5,
                                            L_x_u5);
  vector[n_LON] f_gp_x_u5_LON = PSI_x_u5_LON * (xi_x_u5 .* DELTA_x_u5);
  
  // Term f_gp_x_u6 (GPTerm)
  vector[B_x_u6] DELTA_x_u6 = eigval_shared(alpha_x_u6, ell_x_u6, seq_B_x_u6,
                                            L_x_u6);
  vector[n_LON] f_gp_x_u6_LON = PSI_x_u6_LON * (xi_x_u6 .* DELTA_x_u6);
  
  // Term f_gp_x_u7 (GPTerm)
  vector[B_x_u7] DELTA_x_u7 = eigval_shared(alpha_x_u7, ell_x_u7, seq_B_x_u7,
                                            L_x_u7);
  vector[n_LON] f_gp_x_u7_LON = PSI_x_u7_LON * (xi_x_u7 .* DELTA_x_u7);
  
  // Term f_gp_x_u8 (GPTerm)
  vector[B_x_u8] DELTA_x_u8 = eigval_shared(alpha_x_u8, ell_x_u8, seq_B_x_u8,
                                            L_x_u8);
  vector[n_LON] f_gp_x_u8_LON = PSI_x_u8_LON * (xi_x_u8 .* DELTA_x_u8);
  
  // Term f_gp_x_u9 (GPTerm)
  vector[B_x_u9] DELTA_x_u9 = eigval_shared(alpha_x_u9, ell_x_u9, seq_B_x_u9,
                                            L_x_u9);
  vector[n_LON] f_gp_x_u9_LON = PSI_x_u9_LON * (xi_x_u9 .* DELTA_x_u9);
  
  // Term f_gp_x_u10 (GPTerm)
  vector[B_x_u10] DELTA_x_u10 = eigval_shared(alpha_x_u10, ell_x_u10,
                                              seq_B_x_u10, L_x_u10);
  vector[n_LON] f_gp_x_u10_LON = PSI_x_u10_LON * (xi_x_u10 .* DELTA_x_u10);
  
  // Term f_gp_x_u11 (GPTerm)
  vector[B_x_u11] DELTA_x_u11 = eigval_shared(alpha_x_u11, ell_x_u11,
                                              seq_B_x_u11, L_x_u11);
  vector[n_LON] f_gp_x_u11_LON = PSI_x_u11_LON * (xi_x_u11 .* DELTA_x_u11);
  
  // Term f_gp_x_u12 (GPTerm)
  vector[B_x_u12] DELTA_x_u12 = eigval_shared(alpha_x_u12, ell_x_u12,
                                              seq_B_x_u12, L_x_u12);
  vector[n_LON] f_gp_x_u12_LON = PSI_x_u12_LON * (xi_x_u12 .* DELTA_x_u12);
  
  // Term f_gp_x_u13 (GPTerm)
  vector[B_x_u13] DELTA_x_u13 = eigval_shared(alpha_x_u13, ell_x_u13,
                                              seq_B_x_u13, L_x_u13);
  vector[n_LON] f_gp_x_u13_LON = PSI_x_u13_LON * (xi_x_u13 .* DELTA_x_u13);
  
  // Term f_gp_x_u14 (GPTerm)
  vector[B_x_u14] DELTA_x_u14 = eigval_shared(alpha_x_u14, ell_x_u14,
                                              seq_B_x_u14, L_x_u14);
  vector[n_LON] f_gp_x_u14_LON = PSI_x_u14_LON * (xi_x_u14 .* DELTA_x_u14);
  
  // Term f_gp_x_u15 (GPTerm)
  vector[B_x_u15] DELTA_x_u15 = eigval_shared(alpha_x_u15, ell_x_u15,
                                              seq_B_x_u15, L_x_u15);
  vector[n_LON] f_gp_x_u15_LON = PSI_x_u15_LON * (xi_x_u15 .* DELTA_x_u15);
  
  // Term f_gp_x_u16 (GPTerm)
  vector[B_x_u16] DELTA_x_u16 = eigval_shared(alpha_x_u16, ell_x_u16,
                                              seq_B_x_u16, L_x_u16);
  vector[n_LON] f_gp_x_u16_LON = PSI_x_u16_LON * (xi_x_u16 .* DELTA_x_u16);
  
  // Term f_gp_ageXz (GPTerm)
  vector[B_ageXz * (G_z - 1)] DELTA_ageXz = eigval_grouped(alpha_ageXz,
                                                           ell_ageXz,
                                                           seq_B_ageXz,
                                                           L_ageXz, G_z);
  vector[n_LON] f_gp_ageXz_LON = PSI_ageXz_LON * (xi_ageXz .* DELTA_ageXz);
  
  // Term f_gp_ageXz_u1 (GPTerm)
  vector[B_ageXz_u1 * (G_z_u1 - 1)] DELTA_ageXz_u1 = eigval_grouped(alpha_ageXz_u1,
                                                                    ell_ageXz_u1,
                                                                    seq_B_ageXz_u1,
                                                                    L_ageXz_u1,
                                                                    G_z_u1);
  vector[n_LON] f_gp_ageXz_u1_LON = PSI_ageXz_u1_LON
                                    * (xi_ageXz_u1 .* DELTA_ageXz_u1);
  
  // Term f_gp_ageXz_u2 (GPTerm)
  vector[B_ageXz_u2 * (G_z_u2 - 1)] DELTA_ageXz_u2 = eigval_grouped(alpha_ageXz_u2,
                                                                    ell_ageXz_u2,
                                                                    seq_B_ageXz_u2,
                                                                    L_ageXz_u2,
                                                                    G_z_u2);
  vector[n_LON] f_gp_ageXz_u2_LON = PSI_ageXz_u2_LON
                                    * (xi_ageXz_u2 .* DELTA_ageXz_u2);
  
  // Term f_gp_ageXz_u3 (GPTerm)
  vector[B_ageXz_u3 * (G_z_u3 - 1)] DELTA_ageXz_u3 = eigval_grouped(alpha_ageXz_u3,
                                                                    ell_ageXz_u3,
                                                                    seq_B_ageXz_u3,
                                                                    L_ageXz_u3,
                                                                    G_z_u3);
  vector[n_LON] f_gp_ageXz_u3_LON = PSI_ageXz_u3_LON
                                    * (xi_ageXz_u3 .* DELTA_ageXz_u3);
  
  // Term f_gp_ageXz_u4 (GPTerm)
  vector[B_ageXz_u4 * (G_z_u4 - 1)] DELTA_ageXz_u4 = eigval_grouped(alpha_ageXz_u4,
                                                                    ell_ageXz_u4,
                                                                    seq_B_ageXz_u4,
                                                                    L_ageXz_u4,
                                                                    G_z_u4);
  vector[n_LON] f_gp_ageXz_u4_LON = PSI_ageXz_u4_LON
                                    * (xi_ageXz_u4 .* DELTA_ageXz_u4);
  
  // Term f_gp_ageXz_u5 (GPTerm)
  vector[B_ageXz_u5 * (G_z_u5 - 1)] DELTA_ageXz_u5 = eigval_grouped(alpha_ageXz_u5,
                                                                    ell_ageXz_u5,
                                                                    seq_B_ageXz_u5,
                                                                    L_ageXz_u5,
                                                                    G_z_u5);
  vector[n_LON] f_gp_ageXz_u5_LON = PSI_ageXz_u5_LON
                                    * (xi_ageXz_u5 .* DELTA_ageXz_u5);
  
  // Term f_gp_ageXz_u6 (GPTerm)
  vector[B_ageXz_u6 * (G_z_u6 - 1)] DELTA_ageXz_u6 = eigval_grouped(alpha_ageXz_u6,
                                                                    ell_ageXz_u6,
                                                                    seq_B_ageXz_u6,
                                                                    L_ageXz_u6,
                                                                    G_z_u6);
  vector[n_LON] f_gp_ageXz_u6_LON = PSI_ageXz_u6_LON
                                    * (xi_ageXz_u6 .* DELTA_ageXz_u6);
  
  // Term f_gp_ageXz_u7 (GPTerm)
  vector[B_ageXz_u7 * (G_z_u7 - 1)] DELTA_ageXz_u7 = eigval_grouped(alpha_ageXz_u7,
                                                                    ell_ageXz_u7,
                                                                    seq_B_ageXz_u7,
                                                                    L_ageXz_u7,
                                                                    G_z_u7);
  vector[n_LON] f_gp_ageXz_u7_LON = PSI_ageXz_u7_LON
                                    * (xi_ageXz_u7 .* DELTA_ageXz_u7);
  
  // Term f_gp_ageXz_u8 (GPTerm)
  vector[B_ageXz_u8 * (G_z_u8 - 1)] DELTA_ageXz_u8 = eigval_grouped(alpha_ageXz_u8,
                                                                    ell_ageXz_u8,
                                                                    seq_B_ageXz_u8,
                                                                    L_ageXz_u8,
                                                                    G_z_u8);
  vector[n_LON] f_gp_ageXz_u8_LON = PSI_ageXz_u8_LON
                                    * (xi_ageXz_u8 .* DELTA_ageXz_u8);
  
  // Term f_gp_ageXz_u9 (GPTerm)
  vector[B_ageXz_u9 * (G_z_u9 - 1)] DELTA_ageXz_u9 = eigval_grouped(alpha_ageXz_u9,
                                                                    ell_ageXz_u9,
                                                                    seq_B_ageXz_u9,
                                                                    L_ageXz_u9,
                                                                    G_z_u9);
  vector[n_LON] f_gp_ageXz_u9_LON = PSI_ageXz_u9_LON
                                    * (xi_ageXz_u9 .* DELTA_ageXz_u9);
  
  // Term f_gp_ageXz_u10 (GPTerm)
  vector[B_ageXz_u10 * (G_z_u10 - 1)] DELTA_ageXz_u10 = eigval_grouped(alpha_ageXz_u10,
                                                                    ell_ageXz_u10,
                                                                    seq_B_ageXz_u10,
                                                                    L_ageXz_u10,
                                                                    G_z_u10);
  vector[n_LON] f_gp_ageXz_u10_LON = PSI_ageXz_u10_LON
                                     * (xi_ageXz_u10 .* DELTA_ageXz_u10);
  
  // Term f_gp_ageXz_u11 (GPTerm)
  vector[B_ageXz_u11 * (G_z_u11 - 1)] DELTA_ageXz_u11 = eigval_grouped(alpha_ageXz_u11,
                                                                    ell_ageXz_u11,
                                                                    seq_B_ageXz_u11,
                                                                    L_ageXz_u11,
                                                                    G_z_u11);
  vector[n_LON] f_gp_ageXz_u11_LON = PSI_ageXz_u11_LON
                                     * (xi_ageXz_u11 .* DELTA_ageXz_u11);
  
  // Term f_gp_ageXz_u12 (GPTerm)
  vector[B_ageXz_u12 * (G_z_u12 - 1)] DELTA_ageXz_u12 = eigval_grouped(alpha_ageXz_u12,
                                                                    ell_ageXz_u12,
                                                                    seq_B_ageXz_u12,
                                                                    L_ageXz_u12,
                                                                    G_z_u12);
  vector[n_LON] f_gp_ageXz_u12_LON = PSI_ageXz_u12_LON
                                     * (xi_ageXz_u12 .* DELTA_ageXz_u12);
  
  // Term f_gp_ageXz_u13 (GPTerm)
  vector[B_ageXz_u13 * (G_z_u13 - 1)] DELTA_ageXz_u13 = eigval_grouped(alpha_ageXz_u13,
                                                                    ell_ageXz_u13,
                                                                    seq_B_ageXz_u13,
                                                                    L_ageXz_u13,
                                                                    G_z_u13);
  vector[n_LON] f_gp_ageXz_u13_LON = PSI_ageXz_u13_LON
                                     * (xi_ageXz_u13 .* DELTA_ageXz_u13);
  
  // Term f_gp_ageXz_u14 (GPTerm)
  vector[B_ageXz_u14 * (G_z_u14 - 1)] DELTA_ageXz_u14 = eigval_grouped(alpha_ageXz_u14,
                                                                    ell_ageXz_u14,
                                                                    seq_B_ageXz_u14,
                                                                    L_ageXz_u14,
                                                                    G_z_u14);
  vector[n_LON] f_gp_ageXz_u14_LON = PSI_ageXz_u14_LON
                                     * (xi_ageXz_u14 .* DELTA_ageXz_u14);
  
  // Term f_gp_ageXz_u15 (GPTerm)
  vector[B_ageXz_u15 * (G_z_u15 - 1)] DELTA_ageXz_u15 = eigval_grouped(alpha_ageXz_u15,
                                                                    ell_ageXz_u15,
                                                                    seq_B_ageXz_u15,
                                                                    L_ageXz_u15,
                                                                    G_z_u15);
  vector[n_LON] f_gp_ageXz_u15_LON = PSI_ageXz_u15_LON
                                     * (xi_ageXz_u15 .* DELTA_ageXz_u15);
  
  // Term f_gp_ageXz_u16 (GPTerm)
  vector[B_ageXz_u16 * (G_z_u16 - 1)] DELTA_ageXz_u16 = eigval_grouped(alpha_ageXz_u16,
                                                                    ell_ageXz_u16,
                                                                    seq_B_ageXz_u16,
                                                                    L_ageXz_u16,
                                                                    G_z_u16);
  vector[n_LON] f_gp_ageXz_u16_LON = PSI_ageXz_u16_LON
                                     * (xi_ageXz_u16 .* DELTA_ageXz_u16);
  
  // Sum the terms
  vector[n_LON] f_sum_LON = f_baseline_id_LON + f_gp_age_LON + f_gp_x_LON
                            + f_gp_x_u1_LON + f_gp_x_u2_LON + f_gp_x_u3_LON
                            + f_gp_x_u4_LON + f_gp_x_u5_LON + f_gp_x_u6_LON
                            + f_gp_x_u7_LON + f_gp_x_u8_LON + f_gp_x_u9_LON
                            + f_gp_x_u10_LON + f_gp_x_u11_LON
                            + f_gp_x_u12_LON + f_gp_x_u13_LON
                            + f_gp_x_u14_LON + f_gp_x_u15_LON
                            + f_gp_x_u16_LON + f_gp_ageXz_LON
                            + f_gp_ageXz_u1_LON + f_gp_ageXz_u2_LON
                            + f_gp_ageXz_u3_LON + f_gp_ageXz_u4_LON
                            + f_gp_ageXz_u5_LON + f_gp_ageXz_u6_LON
                            + f_gp_ageXz_u7_LON + f_gp_ageXz_u8_LON
                            + f_gp_ageXz_u9_LON + f_gp_ageXz_u10_LON
                            + f_gp_ageXz_u11_LON + f_gp_ageXz_u12_LON
                            + f_gp_ageXz_u13_LON + f_gp_ageXz_u14_LON
                            + f_gp_ageXz_u15_LON + f_gp_ageXz_u16_LON;
  
  // Gaussian likelihood
  vector[n_LON] h_LON = f_sum_LON + c_hat;
  real log_lik_lon = normal_lpdf(dat_y | h_LON, sigma);
}
model {
  // LonModel
  
  // TermList for f_sum
  // Term f_baseline_id (GroupedOffsetTerm)
  offset_id ~ student_t(4, 0, 5);
  
  // Term f_gp_age (GPTerm)
  alpha_age ~ student_t(20, 0, 1);
  ell_age ~ lognormal(0, 1);
  xi_age ~ normal(0, 1);
  
  // Term f_gp_x (GPTerm)
  alpha_x ~ student_t(20, 0, 1);
  ell_x ~ lognormal(0, 1);
  xi_x ~ normal(0, 1);
  
  // Term f_gp_x_u1 (GPTerm)
  alpha_x_u1 ~ student_t(20, 0, 1);
  ell_x_u1 ~ lognormal(0, 1);
  xi_x_u1 ~ normal(0, 1);
  
  // Term f_gp_x_u2 (GPTerm)
  alpha_x_u2 ~ student_t(20, 0, 1);
  ell_x_u2 ~ lognormal(0, 1);
  xi_x_u2 ~ normal(0, 1);
  
  // Term f_gp_x_u3 (GPTerm)
  alpha_x_u3 ~ student_t(20, 0, 1);
  ell_x_u3 ~ lognormal(0, 1);
  xi_x_u3 ~ normal(0, 1);
  
  // Term f_gp_x_u4 (GPTerm)
  alpha_x_u4 ~ student_t(20, 0, 1);
  ell_x_u4 ~ lognormal(0, 1);
  xi_x_u4 ~ normal(0, 1);
  
  // Term f_gp_x_u5 (GPTerm)
  alpha_x_u5 ~ student_t(20, 0, 1);
  ell_x_u5 ~ lognormal(0, 1);
  xi_x_u5 ~ normal(0, 1);
  
  // Term f_gp_x_u6 (GPTerm)
  alpha_x_u6 ~ student_t(20, 0, 1);
  ell_x_u6 ~ lognormal(0, 1);
  xi_x_u6 ~ normal(0, 1);
  
  // Term f_gp_x_u7 (GPTerm)
  alpha_x_u7 ~ student_t(20, 0, 1);
  ell_x_u7 ~ lognormal(0, 1);
  xi_x_u7 ~ normal(0, 1);
  
  // Term f_gp_x_u8 (GPTerm)
  alpha_x_u8 ~ student_t(20, 0, 1);
  ell_x_u8 ~ lognormal(0, 1);
  xi_x_u8 ~ normal(0, 1);
  
  // Term f_gp_x_u9 (GPTerm)
  alpha_x_u9 ~ student_t(20, 0, 1);
  ell_x_u9 ~ lognormal(0, 1);
  xi_x_u9 ~ normal(0, 1);
  
  // Term f_gp_x_u10 (GPTerm)
  alpha_x_u10 ~ student_t(20, 0, 1);
  ell_x_u10 ~ lognormal(0, 1);
  xi_x_u10 ~ normal(0, 1);
  
  // Term f_gp_x_u11 (GPTerm)
  alpha_x_u11 ~ student_t(20, 0, 1);
  ell_x_u11 ~ lognormal(0, 1);
  xi_x_u11 ~ normal(0, 1);
  
  // Term f_gp_x_u12 (GPTerm)
  alpha_x_u12 ~ student_t(20, 0, 1);
  ell_x_u12 ~ lognormal(0, 1);
  xi_x_u12 ~ normal(0, 1);
  
  // Term f_gp_x_u13 (GPTerm)
  alpha_x_u13 ~ student_t(20, 0, 1);
  ell_x_u13 ~ lognormal(0, 1);
  xi_x_u13 ~ normal(0, 1);
  
  // Term f_gp_x_u14 (GPTerm)
  alpha_x_u14 ~ student_t(20, 0, 1);
  ell_x_u14 ~ lognormal(0, 1);
  xi_x_u14 ~ normal(0, 1);
  
  // Term f_gp_x_u15 (GPTerm)
  alpha_x_u15 ~ student_t(20, 0, 1);
  ell_x_u15 ~ lognormal(0, 1);
  xi_x_u15 ~ normal(0, 1);
  
  // Term f_gp_x_u16 (GPTerm)
  alpha_x_u16 ~ student_t(20, 0, 1);
  ell_x_u16 ~ lognormal(0, 1);
  xi_x_u16 ~ normal(0, 1);
  
  // Term f_gp_ageXz (GPTerm)
  alpha_ageXz ~ student_t(20, 0, 1);
  ell_ageXz ~ lognormal(0, 1);
  xi_ageXz ~ normal(0, 1);
  
  // Term f_gp_ageXz_u1 (GPTerm)
  alpha_ageXz_u1 ~ student_t(20, 0, 1);
  ell_ageXz_u1 ~ lognormal(0, 1);
  xi_ageXz_u1 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u2 (GPTerm)
  alpha_ageXz_u2 ~ student_t(20, 0, 1);
  ell_ageXz_u2 ~ lognormal(0, 1);
  xi_ageXz_u2 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u3 (GPTerm)
  alpha_ageXz_u3 ~ student_t(20, 0, 1);
  ell_ageXz_u3 ~ lognormal(0, 1);
  xi_ageXz_u3 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u4 (GPTerm)
  alpha_ageXz_u4 ~ student_t(20, 0, 1);
  ell_ageXz_u4 ~ lognormal(0, 1);
  xi_ageXz_u4 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u5 (GPTerm)
  alpha_ageXz_u5 ~ student_t(20, 0, 1);
  ell_ageXz_u5 ~ lognormal(0, 1);
  xi_ageXz_u5 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u6 (GPTerm)
  alpha_ageXz_u6 ~ student_t(20, 0, 1);
  ell_ageXz_u6 ~ lognormal(0, 1);
  xi_ageXz_u6 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u7 (GPTerm)
  alpha_ageXz_u7 ~ student_t(20, 0, 1);
  ell_ageXz_u7 ~ lognormal(0, 1);
  xi_ageXz_u7 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u8 (GPTerm)
  alpha_ageXz_u8 ~ student_t(20, 0, 1);
  ell_ageXz_u8 ~ lognormal(0, 1);
  xi_ageXz_u8 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u9 (GPTerm)
  alpha_ageXz_u9 ~ student_t(20, 0, 1);
  ell_ageXz_u9 ~ lognormal(0, 1);
  xi_ageXz_u9 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u10 (GPTerm)
  alpha_ageXz_u10 ~ student_t(20, 0, 1);
  ell_ageXz_u10 ~ lognormal(0, 1);
  xi_ageXz_u10 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u11 (GPTerm)
  alpha_ageXz_u11 ~ student_t(20, 0, 1);
  ell_ageXz_u11 ~ lognormal(0, 1);
  xi_ageXz_u11 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u12 (GPTerm)
  alpha_ageXz_u12 ~ student_t(20, 0, 1);
  ell_ageXz_u12 ~ lognormal(0, 1);
  xi_ageXz_u12 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u13 (GPTerm)
  alpha_ageXz_u13 ~ student_t(20, 0, 1);
  ell_ageXz_u13 ~ lognormal(0, 1);
  xi_ageXz_u13 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u14 (GPTerm)
  alpha_ageXz_u14 ~ student_t(20, 0, 1);
  ell_ageXz_u14 ~ lognormal(0, 1);
  xi_ageXz_u14 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u15 (GPTerm)
  alpha_ageXz_u15 ~ student_t(20, 0, 1);
  ell_ageXz_u15 ~ lognormal(0, 1);
  xi_ageXz_u15 ~ normal(0, 1);
  
  // Term f_gp_ageXz_u16 (GPTerm)
  alpha_ageXz_u16 ~ student_t(20, 0, 1);
  ell_ageXz_u16 ~ lognormal(0, 1);
  xi_ageXz_u16 ~ normal(0, 1);
  
  sigma ~ lognormal(1, 1);
  
  // Likelihood
  if (prior_only == 0) {
    target += log_lik_lon;
  }
}
generated quantities {
  // LonModel
  
  // TermList for f_sum
  // Term f_baseline_id (GroupedOffsetTerm)
  vector[n_LON] f_baseline_id = f_baseline_id_LON;
  
  // Term f_gp_age (GPTerm)
  vector[n_LON] f_gp_age = f_gp_age_LON;
  
  // Term f_gp_x (GPTerm)
  vector[n_LON] f_gp_x = f_gp_x_LON;
  
  // Term f_gp_x_u1 (GPTerm)
  vector[n_LON] f_gp_x_u1 = f_gp_x_u1_LON;
  
  // Term f_gp_x_u2 (GPTerm)
  vector[n_LON] f_gp_x_u2 = f_gp_x_u2_LON;
  
  // Term f_gp_x_u3 (GPTerm)
  vector[n_LON] f_gp_x_u3 = f_gp_x_u3_LON;
  
  // Term f_gp_x_u4 (GPTerm)
  vector[n_LON] f_gp_x_u4 = f_gp_x_u4_LON;
  
  // Term f_gp_x_u5 (GPTerm)
  vector[n_LON] f_gp_x_u5 = f_gp_x_u5_LON;
  
  // Term f_gp_x_u6 (GPTerm)
  vector[n_LON] f_gp_x_u6 = f_gp_x_u6_LON;
  
  // Term f_gp_x_u7 (GPTerm)
  vector[n_LON] f_gp_x_u7 = f_gp_x_u7_LON;
  
  // Term f_gp_x_u8 (GPTerm)
  vector[n_LON] f_gp_x_u8 = f_gp_x_u8_LON;
  
  // Term f_gp_x_u9 (GPTerm)
  vector[n_LON] f_gp_x_u9 = f_gp_x_u9_LON;
  
  // Term f_gp_x_u10 (GPTerm)
  vector[n_LON] f_gp_x_u10 = f_gp_x_u10_LON;
  
  // Term f_gp_x_u11 (GPTerm)
  vector[n_LON] f_gp_x_u11 = f_gp_x_u11_LON;
  
  // Term f_gp_x_u12 (GPTerm)
  vector[n_LON] f_gp_x_u12 = f_gp_x_u12_LON;
  
  // Term f_gp_x_u13 (GPTerm)
  vector[n_LON] f_gp_x_u13 = f_gp_x_u13_LON;
  
  // Term f_gp_x_u14 (GPTerm)
  vector[n_LON] f_gp_x_u14 = f_gp_x_u14_LON;
  
  // Term f_gp_x_u15 (GPTerm)
  vector[n_LON] f_gp_x_u15 = f_gp_x_u15_LON;
  
  // Term f_gp_x_u16 (GPTerm)
  vector[n_LON] f_gp_x_u16 = f_gp_x_u16_LON;
  
  // Term f_gp_ageXz (GPTerm)
  vector[n_LON] f_gp_ageXz = f_gp_ageXz_LON;
  
  // Term f_gp_ageXz_u1 (GPTerm)
  vector[n_LON] f_gp_ageXz_u1 = f_gp_ageXz_u1_LON;
  
  // Term f_gp_ageXz_u2 (GPTerm)
  vector[n_LON] f_gp_ageXz_u2 = f_gp_ageXz_u2_LON;
  
  // Term f_gp_ageXz_u3 (GPTerm)
  vector[n_LON] f_gp_ageXz_u3 = f_gp_ageXz_u3_LON;
  
  // Term f_gp_ageXz_u4 (GPTerm)
  vector[n_LON] f_gp_ageXz_u4 = f_gp_ageXz_u4_LON;
  
  // Term f_gp_ageXz_u5 (GPTerm)
  vector[n_LON] f_gp_ageXz_u5 = f_gp_ageXz_u5_LON;
  
  // Term f_gp_ageXz_u6 (GPTerm)
  vector[n_LON] f_gp_ageXz_u6 = f_gp_ageXz_u6_LON;
  
  // Term f_gp_ageXz_u7 (GPTerm)
  vector[n_LON] f_gp_ageXz_u7 = f_gp_ageXz_u7_LON;
  
  // Term f_gp_ageXz_u8 (GPTerm)
  vector[n_LON] f_gp_ageXz_u8 = f_gp_ageXz_u8_LON;
  
  // Term f_gp_ageXz_u9 (GPTerm)
  vector[n_LON] f_gp_ageXz_u9 = f_gp_ageXz_u9_LON;
  
  // Term f_gp_ageXz_u10 (GPTerm)
  vector[n_LON] f_gp_ageXz_u10 = f_gp_ageXz_u10_LON;
  
  // Term f_gp_ageXz_u11 (GPTerm)
  vector[n_LON] f_gp_ageXz_u11 = f_gp_ageXz_u11_LON;
  
  // Term f_gp_ageXz_u12 (GPTerm)
  vector[n_LON] f_gp_ageXz_u12 = f_gp_ageXz_u12_LON;
  
  // Term f_gp_ageXz_u13 (GPTerm)
  vector[n_LON] f_gp_ageXz_u13 = f_gp_ageXz_u13_LON;
  
  // Term f_gp_ageXz_u14 (GPTerm)
  vector[n_LON] f_gp_ageXz_u14 = f_gp_ageXz_u14_LON;
  
  // Term f_gp_ageXz_u15 (GPTerm)
  vector[n_LON] f_gp_ageXz_u15 = f_gp_ageXz_u15_LON;
  
  // Term f_gp_ageXz_u16 (GPTerm)
  vector[n_LON] f_gp_ageXz_u16 = f_gp_ageXz_u16_LON;
  
  // Sum of terms
  vector[n_LON] f_sum = f_sum_LON;
  
  // Other generated quantities
  vector[n_LON] h = h_LON;
  vector[n_LON] log_lik;
  vector[n_LON] y_log_pred;
  for (i in 1 : n_LON) {
    y_log_pred[i] = normal_rng(h[i], sigma);
    log_lik[i] = normal_lpdf(dat_y[i] | h_LON[i], sigma);
  }
}

