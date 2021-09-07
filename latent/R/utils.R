# REQUIREMENTS ------------------------------------------------------------

# Check if model has correct format
check_model_compatibility <- function(model) {
  stopifnot(is(model, "lgpmodel"))
  ver <- model@info$lgpr_version
  if (ver < "1.1.4" || ver >= "2.0.0") {
    stop(
      "model was created with lgpr ", ver, ", but should be created with ",
      "lgpr version number at least 1.1.4 and less than 2.0.0!"
    )
  }
  TRUE
}

# Check lgpr package version
check_lgpr_version <- function() {
  ver <- packageVersion("lgpr")
  if (ver < "1.1.4" || ver >= "2.0.0") {
    stop(
      "using lgpr version ", ver, ", but ",
      "lgpr version number at least 1.1.4 and less than 2.0.0 is needed!"
    )
  }
  TRUE
}


# MISC --------------------------------------------------------------------

# Get draws of a parameter
get_draws <- function(fit, name) {
  rstan::extract(fit, pars = name)[[name]]
}

# Helper function
compare_param_create_df <- function(fit, fit_approx, name, ag_name) {
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
