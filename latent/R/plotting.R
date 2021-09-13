# Helper function
create_plot_df <- function(data, fit) {
  if (is(fit, "lgpfit")) {
    gpred <- pred(fit)
    f_mean <- as.vector(gpred@f_mean)
    f_sd <- as.vector(gpred@f_std)
  } else if (is(fit, "stanfit")) {
    rv <- as_draws_rvars(fit)$f_latent[1, ]
    f_mean <- as.vector(mean(rv))
    f_sd <- as.vector(sd(rv))
  } else {
    stop("fit should be a stanfit or lgpfit object!")
  }
  cbind(data, f_mean, f_sd)
}

# Plot
plot_f <- function(data, fit, aname = "approx") {
  df <- create_plot_df(data, fit)
  map <- aes(x = age, ymin = f_mean - 2 * f_sd, ymax = f_mean + 2 * f_sd)
  plt <- ggplot(df, aes(x = age, y = f_mean)) +
    geom_ribbon(aes = map, alpha = 0.3) +
    geom_line() +
    ggtitle(aname) +
    ylab("") +
    theme_bw() +
    geom_point(
      data = data, aes(x = age, y = y), pch = 4,
      alpha = 0.3, color = "firebrick3"
    )
  return(plt)
}

# Plot comparison
plot_f_compare_with_exact <- function(data, fit, fit_approx, aname = "approx",
                                      ribbon = FALSE) {
  df <- create_plot_df(data, fit)
  df_approx <- create_plot_df(data, fit_approx)
  N1 <- nrow(df)
  N2 <- nrow(df_approx)
  fit <- as.factor(c(rep("exact", N1), rep(aname, N2)))
  df <- rbind(df, df_approx)
  df <- cbind(df, fit)
  plt <- ggplot(df, aes(x = age, y = f_mean, group = fit, color = fit))
  if (ribbon) {
    plt <- plt + geom_ribbon(
      aes(x = age, ymin = f_mean - 2 * f_sd, ymax = f_mean + 2 * f_sd),
      alpha = 0.15, linetype = 3
    )
  }
  plt <- plt + geom_line() + theme_bw() + ylab("posterior f")
  plt <- plt + geom_point(
    data = data, aes(x = age, y = y), inherit.aes = FALSE,
    pch = 4, alpha = 0.3
  )
  return(plt)
}

# Plot means comparison in same panel
plot_f_compare_same <- function(data, fits) {
  nams <- names(fits)
  J <- length(fits)
  DF <- NULL
  TYP <- NULL
  for (j in 1:J) {
    df <- create_plot_df(data, fits[[j]])
    DF <- rbind(DF, df)
    TYP <- c(TYP, c(rep(nams[j], nrow(df))))
  }
  fit <- as.factor(TYP)
  plt <- ggplot(DF, aes(x = age, y = f_mean, group = fit, color = fit))
  plt <- plt + geom_line()
  plt <- plt + geom_point(
    data = data, aes(x = age, y = y), inherit.aes = FALSE,
    pch = 4, alpha = 0.3
  )
  plt <- plt + theme_bw() + ylab("posterior mean f")
  return(plt)
}

# Plot mean and sd comparison in separate figures
plot_f_compare_separate <- function(data, fits, ncol = 2, nrow = NULL,
                                    last_is_exact = FALSE) {
  PLOTS <- list()
  J <- length(fits)
  nams <- names(fits)
  if (last_is_exact) {
    fit_exact <- fits[[J]]
    fits <- fits[1:(J - 1)]
    J <- length(fits)
  } else {
    fit_exact <- NULL
  }
  for (j in seq_len(J)) {
    if (!(is.null(fit_exact))) {
      plt <- plot_f_compare_with_exact(data, fit_exact, fits[[j]],
        ribbon = T, aname = nams[j]
      )
    } else {
      plt <- plot_f(data, fits[[j]], aname = nams[j])
    }
    PLOTS[[j]] <- plt
  }
  if (is.null(nrow)) nrow <- ceiling(J / ncol)
  ggarrange(plotlist = PLOTS, nrow = nrow, ncol = ncol)
}

# Pairs plot of alpha and ell parameters
pairs_kernelparams <- function(fit) {
  stopifnot(is(fit, "stanfit"))
  pairs(fit[[1]], pars = c("alpha", "ell"))
}

# Compare approximate and exact EQ covariance functions
plot_kernelcomparison_eq <- function(pars_approx, pars, stan_data, idx_x = 1) {
  cmp <- compare_kernels_eq(pars_approx, pars, stan_data, idx_x)
  plot(cmp$x, cmp$K[1, ], type = "l", ylim = c(0, alpha^2), xlab = "x", ylab = "K")
  lines(cmp$x, cmp$K_approx[1, ], col = "firebrick")
}
