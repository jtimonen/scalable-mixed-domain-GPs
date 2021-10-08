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
    fd <- get_cmdstanfit(fit)
    s_draws <- as.vector(posterior::merge_chains(fd$draws("sigma")))
    # scale std to original data scale
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

# Root mean squared error (f marginalized)
compute_rmse.marginal <- function(pred, y_star) {
  stopifnot(is(pred, "GaussianPrediction"))
  mu <- colMeans(pred@y_mean) # mean estimates
  se <- (mu - y_star)**2
  return(sqrt(mean(se)))
}

# Root mean squared error (f sampled)
compute_rmse.sampled <- function(pred, y_star) {
  stopifnot(is(pred, "Prediction"))
  mu <- colMeans(pred@h) # mean estimates
  se <- (mu - y_star)**2
  return(sqrt(mean(se)))
}

# Compute root mean squared error
compute_rmse <- function(fit, pred, y_star) {
  if (isa(fit, "lgpfit")) {
    rmse <- compute_rmse.marginal(pred, y_star)
  } else {
    rmse <- compute_rmse.sampled(pred, y_star)
  }
  return(rmse)
}


# Formatting
format_result_row <- function(res, scales, num_bfs) {
  N1 <- length(scales)
  N2 <- length(num_bfs)
  M <- matrix(res[1:(N1 * N2)], N1, N2, byrow = TRUE)
  rownames(M) <- paste0("c = ", scales)
  colnames(M) <- paste0("B = ", num_bfs)
  list(
    approx = M,
    exact = res[N1 * N2 + 1]
  )
}

# Formatting
format_results <- function(res, scales, num_bfs) {
  out <- list()
  confs <- colnames(res)
  J <- nrow(res)
  for (j in 1:J) {
    out[[j]] <- format_result_row(res[j, ], scales, num_bfs)
  }
  names(out) <- rownames(res)
  out[[J + 1]] <- format_result_row(confs, scales, num_bfs)
  names(out)[J + 1] <- "conf"
  return(out)
}

# Compute all evaluation metrics
compute_metrics <- function(fits, preds, y_star) {
  num_fits <- length(fits)
  elpds_w1 <- rep(0.0, num_fits)
  elpds_w2 <- rep(0.0, num_fits)
  rmses <- rep(0.0, num_fits)
  for (j in 1:num_fits) {
    elpds_w1[j] <- compute_elpd(fits[[j]], preds[[j]], y_star, 1)
    elpds_w2[j] <- compute_elpd(fits[[j]], preds[[j]], y_star, 2)
    rmses[j] <- compute_rmse(fits[[j]], preds[[j]], y_star)
  }
  res <- rbind(elpds_w1, elpds_w2, rmses)
  colnames(res) <- names(fits)
  return(res)
}
