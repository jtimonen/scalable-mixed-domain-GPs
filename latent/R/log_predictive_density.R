# First way
lpd_m_way1 <- function(pred, y_star) {
  S <- lgpr::num_paramsets(pred)
  y_means <- pred@y_mean
  y_stds <- pred@y_std
  log_pds <- array(0.0, dim = dim(y_means))
  for (s in seq_len(S)) {
    mu <- y_means[s, ]
    sig <- y_stds[s, ]
    log_pds[s, ] <- stats::dnorm(y_star, mean = mu, sd = sig, log = TRUE)
  }
  return(rowMeans(log_pds))
}

# Second way
lpd_m_way2 <- function(pred, y_star) {
  mu <- colMeans(pred@y_mean) # mean estimates
  sig2 <- colMeans(pred@y_std**2) # variance estimates
  sig <- sqrt(sig2)
  log_pds <- stats::dnorm(y_star, mean = mu, sd = sig, log = TRUE)
  return(log_pds)
}

# Log predictive density (f marginalized)
compute_lpd.marginal <- function(pred, y_star) {
  stopifnot(is(pred, "GaussianPrediction"))
  list(
    way1 = lpd_m_way1(pred, y_star),
    way2 = lpd_m_way2(pred, y_star)
  )
}


# First way
lpd_sg_way1 <- function(pred, y_star, s_draws) {
  S <- lgpr::num_paramsets(pred)
  log_pds <- array(0.0, dim = dim(pred@h))
  for (s in seq_len(S)) {
    mu <- pred@h[s, ]
    sig <- s_draws[s]
    log_pds[s, ] <- stats::dnorm(y_star, mean = mu, sd = sig, log = TRUE)
  }
  return(rowMeans(log_pds))
}

# Second way
lpd_sg_way2 <- function(pred, y_star, s_draws) {
  S <- lgpr::num_paramsets(pred)
  log_pds <- array(0.0, dim = dim(pred@h))
  sig2 <- mean(s_draws**2) # estimate for sigma2
  sig <- sqrt(sig2)
  h_means <- colMeans(pred@h)
  h_std <- apply(pred@h, 2, stats::sd)
  log_pds <- stats::dnorm(y_star, mean = h_means, sd = h_std + sig, log = TRUE)
  return(log_pds)
}

# Log predictive density (f sampled)
compute_lpd.sampled_gaussian <- function(pred, y_star, s_draws) {
  stopifnot(is(pred, "Prediction"))
  list(
    way1 = lpd_sg_way1(pred, y_star, s_draws),
    way2 = lpd_sg_way2(pred, y_star, s_draws)
  )
}

# Compute expected log predictive density
compute_elpd <- function(fit, pred, y_star, way = 1) {
  if (isa(fit, "lgpfit")) {
    lpd <- compute_lpd.marginal(pred, y_star)
  } else {
    fd <- fit@fit[[1]]
    s_draws <- as.vector(posterior::merge_chains(fd$draws("sigma")))
    # scale variance to original data scale
    y_scl <- fit@model@exact_model@var_scalings$y
    s_draws <- s_draws * y_scl@scale
    lpd <- compute_lpd.sampled_gaussian(pred, y_star, s_draws)
  }
  if (way == 2) {
    field <- "way2"
  } else {
    field <- "way1"
  }
  mean(lpd[[field]])
}
