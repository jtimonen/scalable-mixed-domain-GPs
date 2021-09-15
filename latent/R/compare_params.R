
# Pairs plot of alpha and ell parameters
pairs_kernelparams <- function(fit) {
  stopifnot(is(fit, "stanfit"))
  pairs(fit, pars = c("alpha", "ell"))
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
