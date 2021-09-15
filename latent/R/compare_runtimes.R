
# Get chain runtimes
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
