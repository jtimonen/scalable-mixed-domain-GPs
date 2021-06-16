# Check lgpr version
check_version <- function(pkg_name) {
  cur_ver <- packageVersion("lgpr")
  min_ver <- "1.1"
  if (cur_ver < min_ver) {
    msg <- sprintf(
      "lgpr version at least %s is needed, found %s", min_ver, cur_ver
    )
    stop(msg)
  }
}

# Get draws of a parameter
get_draws <- function(fit, name) {
  df <- as.data.frame(fit$draws(variables=name))
  as.vector(as.matrix(df))
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
compare_param <- function(fit, fit_approx, name, ag_name="approx") {
  df <- compare_param_create_df(fit, fit_approx, name, ag_name)
  mu <- compare_param_means(fit, fit_approx, name, ag_name)
  plt <- ggplot(df, aes(x=value, color=label)) + geom_density() +
    geom_vline(data=mu, aes(xintercept=mean, color=label),
               linetype="dashed")
  plt <- plt + ggtitle(name) + 
    theme(legend.position="top") + theme(legend.title = element_blank())
  return(plt)
}

# Compare runtimes of chains
compare_runtime <- function(fit, fit_approx, ag_name="approx") {
  t <- fit$time()$chains$total
  t_approx <- fit_approx$time()$chains$total
  g <- rep("exact", length(t))
  g_approx <- rep(ag_name, length(t_approx))
  runtime <- c(t, t_approx)
  label <- as.factor(c(g, g_approx))
  data.frame(runtime, label)
}

# Plot runtimes
compare_runtime_plot <- function(fit, fit_approx, ag_name) {
  df <- compare_runtime(fit, fit_approx, ag_name=ag_name)
  plt <- ggplot(df, aes(y=runtime, x=label, fill=label)) + geom_boxplot()
  plt <- plt + ggtitle("Chain runtimes (s)") + 
    theme(legend.title = element_blank())
  return(plt)
}

# Parameter comparison plot
plot_params_comparison <- function(fit, fit_approx, ag_name = "approx") {
  a1 <- compare_param(fit, fit_approx, "alpha[1]", ag_name)
  a2 <- compare_param(fit, fit_approx, "alpha[2]", ag_name)
  e1 <- compare_param(fit, fit_approx, "ell[1]", ag_name)
  e2 <- compare_param(fit, fit_approx, "ell[2]", ag_name)
  s1 <- compare_param(fit, fit_approx, "sigma[1]", ag_name)
  r <- compare_runtime_plot(fit, fit_approx, ag_name=ag_name)
  plots <- ggarrange(a1, a2, e1, e2, s1, r, labels = "auto")
  return(plots)
}
