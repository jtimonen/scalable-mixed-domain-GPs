# Wraps STAN_create_psi_mats
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

# Wraps computing all the transformed data
do_transformed_data <- function(stan_data) {
  expose_stanfuns()
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

# Helper functions
get_runtimes <- function(x) {
  if (is(x, "lgpfit")) x <- x@stan_fit
  if (is(x, "ApproxModelFit")) {
    fit <- get_cmdstanfit(x)
    tims <- fit$time()$chains$total
  } else {
    tims <- as.numeric(rowSums(get_elapsed_time(x)))
  }
  return(tims)
}
get_ndiv <- function(x) {
  if (is(x, "lgpfit")) x <- x@stan_fit
  if (is(x, "ApproxModelFit")) {
    fit <- x@fit[[1]]
    ndiv <- sum(fit$sampler_diagnostics()[, , "divergent__"])
  } else {
    ndiv <- sum(rstan::get_divergent_iterations(x))
  }
  return(ndiv)
}

# Get param draws
get_pars <- function(x) {
  pars <- c("alpha", "ell", "sigma")
  if (is(x, "lgpfit")) x <- x@stan_fit
  if (is(x, "ApproxModelFit")) {
    fit <- x@fit[[1]]
    d <- fit$draws(pars)
  } else {
    d <- posterior::subset_draws(posterior::as_draws(x), pars)
  }
  d <- posterior::merge_chains(d)
  return(d)
}

# Get experiment results
summarize_results <- function(fits) {
  draws <- lapply(fits, get_pars)
  p_means <- sapply(draws, function(x) apply(x, 3, mean))
  p_sds <- sapply(draws, function(x) apply(x, 3, stats::sd))
  colnames(p_means) <- names(fits)
  colnames(p_sds) <- names(fits)
  list(
    p_means = p_means,
    p_sds = p_sds,
    runtimes = sapply(fits, get_runtimes),
    num_div = sapply(fits, get_ndiv)
  )
}
