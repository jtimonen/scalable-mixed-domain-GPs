# Log predictive density (f marginalized)
compute_lpd.marginal <- function(pred, y_star) {
  stopifnot(is(pred, "GaussianPrediction"))
  P <- lgpr::num_evalpoints(pred)
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
  # sig <- mean(s_draws) # estimate for sigma
  # mu <- colMeans(pred@f) # estimate for f
  # log_pds <- stats::dnorm(y_star, mean = mu, sd = sig, log = TRUE)
}

# Log predictive density (f sampled)
compute_lpd.sampled_gaussian <- function(pred, y_star, s_draws) {
  stopifnot(is(pred, "Prediction"))
  P <- lgpr::num_evalpoints(pred)
  S <- lgpr::num_paramsets(pred)
  y_means <- pred@f
  log_pds <- array(0.0, dim = dim(y_means))
  sig <- mean(s_draws) # estimate for sigma
  h_means <- colMeans(pred@h)
  h_std <- apply(pred@h, 2, stats::sd)
  # for (s in seq_len(S)) {
  #  mu <- y_means[s, ]
  #  sig <- s_draws[s]
  #  log_pds[s, ] <- stats::dnorm(y_star, mean = mu, sd = sig, log = TRUE)
  # }
  # return(rowMeans(log_pds))
  log_pds <- stats::dnorm(y_star, mean = h_means, sd = h_std + sig, log = TRUE)
  return(log_pds)
}


# Compute expected log predictive density
compute_elpd <- function(fit, pred, y_star) {
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
  mean(lpd)
}
