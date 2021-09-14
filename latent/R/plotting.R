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
  pairs(fit, pars = c("alpha", "ell"))
}

# Compare approximate and exact EQ covariance functions
plot_kernelcomparison_eq <- function(pars_approx, pars, stan_data, idx_x = 1) {
  cmp <- compare_kernels_eq(pars_approx, pars, stan_data, idx_x)
  alpha <- max(pars_approx[1], pars[1])
  r <- cmp$x - cmp$x[1]
  plot(r, cmp$K[1, ],
    type = "l", ylim = c(0, alpha^2),
    xlab = "r", ylab = "K", main = "Covariance function (EQ)"
  )
  lines(r, cmp$K_approx[1, ], col = "firebrick")
}


# Create a data frame for ggplot
compare_param_create_df <- function(fit, fit_approx, name, ag_name) {
  get_draws <- function(fit, name) {
    rstan::extract(fit, pars = name)[[name]]
  }
  d <- get_draws(fit, name)
  d_approx <- get_draws(fit_approx, name)
  g <- rep("exact", length(d))
  g_approx <- rep(ag_name, length(d_approx))
  value <- c(d, d_approx)
  label <- as.factor(c(g, g_approx))
  data.frame(value, label)
}

# Helper function
compare_param_means <- function(fit, fit_approx, name, ag_name) {
  d <- get_draws(fit, name)
  d_approx <- get_draws(fit_approx, name)
  m <- mean(d)
  m_approx <- mean(d_approx)
  mean <- c(m, m_approx)
  label <- c("exact", ag_name)
  data.frame(mean, label)
}

# Compare draws from exact and approximate fit
compare_param <- function(fit, fit_approx, name, ag_name = "approx") {
  df <- compare_param_create_df(fit, fit_approx, name, ag_name)
  mu <- compare_param_means(fit, fit_approx, name, ag_name)
  plt <- ggplot(df, aes(x = value, color = label)) +
    geom_density() +
    geom_vline(
      data = mu, aes(xintercept = mean, color = label),
      linetype = "dashed"
    )
  plt <- plt + ggtitle(name) +
    theme(legend.position = "top") + theme(legend.title = element_blank())
  return(plt)
}

get_chain_times <- function(fit) {
  as.vector(rowSums(get_elapsed_time(fit)))
}

# Compare runtimes of chains
compare_runtime <- function(fit, fit_approx, ag_name = "approx") {
  t <- get_chain_times(fit)
  t_approx <- get_chain_times(fit_approx)
  m <- round(mean(t), 2)
  m_approx <- round(mean(t_approx), 1)
  RED <- "\u001b[31m"
  GREEN <- "\u001b[32m"
  if (m_approx > m) {
    C_approx <- RED
    C_exact <- GREEN
  } else {
    C_approx <- GREEN
    C_exact <- RED
  }
  NC <- "\u001b[0m"
  s1 <- paste0(C_exact, m, " (exact)", NC)
  s2 <- paste0(C_approx, m_approx, " (approx)", NC)
  cat("Average chain runtimes (s): ", s1, " vs. ", s2, "\n", sep = "")
  g <- rep("exact", length(t))
  g_approx <- rep(ag_name, length(t_approx))
  runtime <- c(t, t_approx)
  label <- as.factor(c(g, g_approx))
  data.frame(runtime, label)
}

# Plot runtimes
compare_runtime_plot <- function(fit, fit_approx, ag_name, N) {
  df <- compare_runtime(fit, fit_approx, ag_name = ag_name)
  plt <- ggplot(df, aes(y = runtime, x = label)) +
    geom_boxplot()
  plt <- plt + theme(legend.title = element_blank())
  ymax <- max(df$runtime)
  plt <- plt + ggtitle(paste0("num_obs = ", N))
  plt <- plt + coord_cartesian(ylim = c(0, ymax + 10)) + ylab("runtime (s)")
  return(plt)
}

# Parameter comparison plot
plot_params_comparison <- function(fit, fit_approx, ag_name = "approx", N = 0) {
  a1 <- compare_param(fit, fit_approx, "alpha[1]", ag_name)
  a2 <- compare_param(fit, fit_approx, "alpha[2]", ag_name)
  a3 <- compare_param(fit, fit_approx, "alpha[3]", ag_name)
  e1 <- compare_param(fit, fit_approx, "ell[1]", ag_name)
  e2 <- compare_param(fit, fit_approx, "ell[2]", ag_name)
  e3 <- compare_param(fit, fit_approx, "ell[3]", ag_name)
  s1 <- compare_param(fit, fit_approx, "sigma[1]", ag_name)
  r <- compare_runtime_plot(fit, fit_approx, ag_name = ag_name, N = N)
  plots <- ggarrange(a1, a2, a3, e1, e2, e3, s1, r, labels = "auto")
  return(plots)
}

# Parameter comparison plot
plot_params_comparison_onecomp <- function(fit, fit_approx, ag_name = "approx",
                                           N = 0) {
  a1 <- compare_param(fit, fit_approx, "alpha[1]", ag_name)
  e1 <- compare_param(fit, fit_approx, "ell[1]", ag_name)
  s1 <- compare_param(fit, fit_approx, "sigma[1]", ag_name)
  r <- compare_runtime_plot(fit, fit_approx, ag_name = ag_name, N = N)
  plots <- ggarrange(a1, e1, s1, r, labels = "auto")
  return(plots)
}
