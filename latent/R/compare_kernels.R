# Approximate the EQ kernel
approximate_kernel_eq <- function(alpha, ell, stan_data, idx_x = 1) {
  td <- do_transformed_data(stan_data) # exposes stan functions
  PHI <- td$PHI_mats[[idx_x]]
  dj <- STAN_basisfun_eq_multipliers(alpha, ell, td$seq_B, td$L[idx_x])
  DELTA <- diag(dj)
  K <- PHI %*% DELTA %*% t(PHI)
  return(K)
}

# Compare approximate and exact EQ covariance functions
compare_kernels_eq <- function(pars, stan_dats, idx_x = 1, exact_name) {
  J <- ncol(pars)
  x <- stan_dats[[1]]$x_cont[idx_x, ]
  r <- x - x[1]
  R <- c()
  K <- c()
  model <- c()
  for (j in 1:J) {
    nam <- colnames(pars)[j]
    pars_j <- pars[, j]
    alpha <- as.numeric(pars_j["alpha[1]"])
    ell <- as.numeric(pars_j["ell[1]"])
    if (nam == exact_name) {
      k <- STAN_kernel_eq(x, x, alpha, ell)
    } else {
      k <- approximate_kernel_eq(alpha, ell, stan_dats[[j]], idx_x)
    }
    R <- c(R, r)
    K <- c(K, k[1, ])
    model <- c(model, rep(nam, length(r)))
  }
  df <- data.frame(R, K, as.factor(model))
  colnames(df) <- c("r", "k", "model")
  return(df)
}

# Compare approximate and exact EQ covariance functions
plot_kernelcomparison_eq <- function(pars, stan_data, idx_x = 1,
                                     exact_name = "marginal") {
  df <- compare_kernels_eq(pars, stan_data, idx_x, exact_name)
  plt <- ggplot(df, aes(x = r, y = k, group = model, color = model)) +
    geom_line()
  plt <- plt +
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    )
  return(plt)
}
