# Matrix rows to a list
#
# @param x a matrix or array with \code{m} rows and \code{n} columns
# @return an unnamed list of length \code{m} where each element is a
# vector of length \code{n}
matrix_to_list <- function(x) {
  m <- dim(x)[1]
  L <- list()
  for (i in seq_len(m)) {
    L[[i]] <- x[i, ]
  }
  return(L)
}


# Call STAN_create_phi_mats
create_phi_mats <- function(stan_data) {
  N <- stan_data$num_obs
  D <- stan_data$num_cov_cont
  seq_M <- STAN_seq_len(stan_data$num_bf)
  c <- stan_data$scale_bf
  x_cont <- matrix_to_list(stan_data$x_cont)
  X_hr <- as.array(stan_data$X_hr)
  PHI <- STAN_create_phi_mats(N, D, seq_M, c, x_cont, X_hr)
  names(PHI) <- rownames(x_cont)
  return(PHI)
}

# Call STAN_create_psi_mats
create_psi_mats <- function(stan_data, PHI_mats) {
  N <- stan_data$num_obs
  M <- stan_data$num_bf
  num_xi <- as.array(stan_data$num_xi)
  comps <- matrix_to_list(stan_data$components)
  x_cat <- matrix_to_list(stan_data$x_cat)
  C_ranks <- as.array(stan_data$C_ranks)
  C_sizes <- as.array(stan_data$C_sizes)
  C_rsp <- as.array(stan_data$C_rsp)
  C_vals <- as.array(stan_data$C_vals)
  C_vecs <- as.array(stan_data$C_vecs)
  STAN_create_psi_mats(
    N, M, num_xi, PHI_mats, comps,
    x_cat, C_vals, C_vecs, C_ranks, C_sizes, C_rsp
  )
}
