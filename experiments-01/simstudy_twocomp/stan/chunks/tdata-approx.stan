  vector[num_bf] seq_B = STAN_seq_len(num_bf);
  matrix[num_obs, num_bf] mat_B = 
    transpose(rep_matrix(seq_B, num_obs));
  vector[num_cov_cont] L = scale_bf * X_hr;
  matrix[num_obs, num_bf] PHI_mats[num_cov_cont] = 
    STAN_create_basisfun_mats(x_cont, mat_B, L);
  matrix[num_obs, sum(num_xi)] PSI_mats = STAN_create_psi_mats(num_xi, PHI_mats,
    components, x_cat, C_vals, C_vecs, C_ranks, C_sizes, C_rsp);
