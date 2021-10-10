# Log predictive density (using f samples)
compute_lpd.sampled_gaussian <- function(pred, y_star, s_draws) {
  stopifnot(is(pred, "Prediction"))
  S <- lgpr::num_paramsets(pred)
  log_pds <- array(0.0, dim = dim(pred@h))
  for (s in seq_len(S)) {
    mu <- pred@h[s, ]
    sig <- s_draws[s]
    log_pds[s, ] <- stats::dnorm(y_star, mean = mu, sd = sig, log = TRUE)
  }
  return(rowMeans(log_pds))
}

# Get draws of sigma
get_sigma_draws <- function(fit) {
  if (is(fit, "ApproxModelFit")) {
    fd <- get_cmdstanfit(fit)
    s_draws <- posterior::merge_chains(fd$draws("sigma"))
  } else {
    s_draws <- lgpr::get_draws(fit, pars = "sigma")
  }
  return(as.vector(s_draws))
}

# Compute mean log predictive density
compute_mlpd <- function(fit, pred, y_star, way = 1) {
  stopifnot(is(pred, "Prediction"))
  s_draws <- get_sigma_draws(fit)
  # scale std to original data scale
  if (is(fit, "lgpfit")) {
    emodel <- fit@model
  } else {
    emodel <- fit@model@exact_model
  }
  y_scl <- emodel@var_scalings$y
  s_draws <- s_draws * y_scl@scale
  lpd <- compute_lpd.sampled_gaussian(pred, y_star, s_draws)
  mean(lpd)
}

# Compute root mean squared error
compute_rmse <- function(fit, pred, y_star) {
  stopifnot(is(pred, "Prediction"))
  mu <- colMeans(pred@h) # mean estimates
  se <- (mu - y_star)**2
  return(sqrt(mean(se)))
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
  mlpd <- rep(0.0, num_fits)
  rmse <- rep(0.0, num_fits)
  for (j in 1:num_fits) {
    mlpd[j] <- compute_mlpd(fits[[j]], preds[[j]], y_star)
    rmse[j] <- compute_rmse(fits[[j]], preds[[j]], y_star)
  }
  res <- rbind(mlpd, rmse)
  colnames(res) <- names(fits)
  return(res)
}
