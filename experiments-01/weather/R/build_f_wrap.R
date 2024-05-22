# Wraps STAN_build_f_latent
build_f <- function(stan_data, tdata, alpha, ell, xi) {
  expose_stanfuns()
  N <- stan_data$num_obs
  seq_M <- STAN_seq_len(stan_data$num_bf)
  scale_bf <- stan_data$scale_bf
  num_xi <- as.array(stan_data$num_xi)
  comps <- matrix_to_list(stan_data$components)
  C_ranks <- as.array(stan_data$C_ranks)
  STAN_build_f(
    comps, num_xi, C_ranks, seq_M, tdata$L, tdata$PSI, alpha, ell, xi
  )
}

# Wraps STAN_build_f_draws
build_f_draws <- function(stan_data, tdata, ALPHA, ELL, XI) {
  expose_stanfuns()
  N <- stan_data$num_obs
  seq_M <- STAN_seq_len(stan_data$num_bf)
  scale_bf <- stan_data$scale_bf
  num_xi <- as.array(stan_data$num_xi)
  comps <- matrix_to_list(stan_data$components)
  C_ranks <- as.array(stan_data$C_ranks)
  STAN_build_f_draws(
    comps, num_xi, C_ranks, seq_M, tdata$L, tdata$PSI, ALPHA, ELL, XI
  )
}
