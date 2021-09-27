# Compute log predictive density
compute_lpd <- function(fit, df_star) {
  stopifnot(is(fit, "lgpfit"))
  p <- lgpr::pred(fit, x = df_star, reduce = NULL)
  y_name <- fit@model@var_names$y
  y_star <- df_star[[y_name]]
  gaussian_lpd(p, y_star)
}

# Gaussian log predictive density
gaussian_lpd <- function(pred, y_star) {
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
  return(rowSums(log_pds))
}
