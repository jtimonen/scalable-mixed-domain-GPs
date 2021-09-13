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
    num_xi, PHI_mats, comps,
    x_cat, C_vals, C_vecs, C_ranks, C_sizes, C_rsp
  )
}

# Compute all the transformed data in R, by calling exposed Stan functions
do_transformed_data <- function(stan_data) {
  num_bf <- stan_data$num_bf
  num_obs <- stan_data$num_obs
  seq_B <- STAN_seq_len(num_bf)
  X_hr <- as.vector(stan_data$X_hr)
  mat_B <- matrix(rep(seq_B, each = num_obs), num_obs, num_bf, byrow = FALSE)
  L <- stan_data$scale_bf * X_hr
  x_cont <- matrix_to_list(stan_data$x_cont)
  PHI_mats <- STAN_create_basisfun_mats(x_cont, mat_B, L)
  PSI_mats <- create_psi_mats(stan_data, PHI_mats)
  list(
    seq_B = seq_B,
    mat_B = mat_B,
    L = L,
    PHI_mats = PHI_mats,
    PSI_mats = PSI_mats
  )
}

# Call STAN_build_f_latent
build_f_latent <- function(stan_data, PSI, alpha, ell, xi) {
  N <- stan_data$num_obs
  seq_M <- STAN_seq_len(stan_data$num_bf)
  scale_bf <- stan_data$scale_bf
  X_hr <- as.array(stan_data$X_hr)
  num_xi <- as.array(stan_data$num_xi)
  comps <- matrix_to_list(stan_data$components)
  C_ranks <- as.array(stan_data$C_ranks)
  STAN_build_f_latent(
    N, comps, num_xi, C_ranks, seq_M, X_hr, scale_bf, PSI,
    alpha, ell, xi
  )
}

# Draw xi and build latent
draw_f_latent <- function(stan_data, PSI, alpha = NULL, ell = NULL) {
  if (is.null(alpha)) {
    alpha <- rep(1.0, stan_data$num_comps)
  }
  if (is.null(ell)) {
    ell <- rep(1.0, stan_data$num_ell)
  }
  alpha <- as.array(alpha)
  ell <- as.array(ell)
  P <- sum(stan_data$num_xi)
  xi <- rnorm(n = P)
  build_f_latent(stan_data, PSI, alpha, ell, xi)
}

# Approximate the EQ kernel
approximate_kernel_eq <- function(alpha, ell, stan_data, idx_x = 1) {
  td <- do_transformed_data(stan_data)
  PHI <- td$PHI_mats[[idx_x]]
  dj <- STAN_basisfun_eq_multipliers(alpha, ell, td$seq_B, td$L[idx_x])
  DELTA <- diag(dj)
  K <- PHI %*% DELTA %*% t(PHI)
  list(K_approx = K, x = stan_data$x_cont[idx_x, ])
}
