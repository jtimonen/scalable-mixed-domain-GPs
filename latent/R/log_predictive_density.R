# Compute log predictive density
compute_lpd <- function(fit, df_star) {
  stopifnot(is(fit, "lgpfit"))
  p <- lgpr::pred(fit, x = df_star, reduce = NULL)
  y_name <- fit@model@var_names$y
  y_star <- df_star[[y_name]]
  gaussian_lpd(p, y_star)
}

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
}

# Log predictive density (f sampled)
compute_lpd.sampled_gaussian <- function(pred, y_star, s_draws) {
  stopifnot(is(pred, "Prediction"))
  P <- lgpr::num_evalpoints(pred)
  S <- lgpr::num_paramsets(pred)
  y_means <- pred@f
  log_pds <- array(0.0, dim = dim(y_means))
  for (s in seq_len(S)) {
    mu <- y_means[s, ]
    sig <- s_draws[s]
    log_pds[s, ] <- stats::dnorm(y_star, mean = mu, sd = sig, log = TRUE)
  }
  return(rowMeans(log_pds))
}


# Compute expected log predictive density
compute_elpd <- function(model, fit, df_star, num_bf = NULL, scale_bf = NULL) {
  y_name <- lgpr:::get_y_name(model)
  y_star <- df_star[[y_name]]
  if (isa(fit, "lgpfit")) {
    p <- lgpr::pred(fit, x = df_star, reduce = NULL)
    lpd <- compute_lpd.marginal(p, y_star)
  } else {
    p_approx <- pred_approx(model, fit, df_star, num_bf, scale_bf)
    s_draws <- as.vector(posterior::merge_chains(fit$draws("sigma")))
    lpd <- compute_lpd.sampled_gaussian(p_approx, y_star, s_draws)
  }
  return(mean(lpd))
}
