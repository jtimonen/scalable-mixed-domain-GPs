
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

# Compare runtimes
plot_runtimes <- function(PRES, NUM_BF, N) {
  tims <- sapply(PRES, function(x) {
    getElement(x, "t_means")
  })
  t_exact <- mean(tims["marginal", ])
  tims <- tims[1:nrow(tims) - 1, ]
  cc <- as.factor(rep(colnames(tims), each = nrow(tims)))
  bb <- rep(NUM_BF, times = ncol(tims))
  df <- data.frame(as.vector(tims), cc, bb)
  main <- paste0("N = ", N)
  sub <- paste0("average chain time for exact fit = ", t_exact, " s")
  colnames(df) <- c("time", "c", "B")
  plt <- ggplot(df, aes(x = bb, y = time, group = c, color = c)) +
    geom_line() +
    geom_point()
  plt <- plt + ggtitle(main, subtitle = sub) + ylab(" average chain time (s)") +
    xlab("number of basis functions")
  return(plt)
}


# Compare runtimes
plot_runtimes_wrt_N <- function(PRES, NUM_BF, N_sizes, scale_bf) {
  tims <- sapply(PRES, function(x) {
    getElement(x, "t_means")
  })
  t_exact <- tims["marginal", ]
  tims <- tims[1:nrow(tims) - 1, ]
  nn <- rep(N_sizes, each = nrow(tims))
  bb <- as.factor(rep(NUM_BF, times = ncol(tims)))
  df <- data.frame(as.vector(tims), nn, bb)
  main <- paste0("c = ", scale_bf)
  sub <- paste0("black dotted  line is for exact fit")
  colnames(df) <- c("time", "N", "B")
  plt <- ggplot(df, aes(x = N, y = time, group = B, color = B)) +
    geom_line() +
    geom_point()
  plt <- plt + ggtitle(main, subtitle = sub) + ylab(" average chain time (s)") +
    xlab("N")
  df_exact <- data.frame(N=N_sizes[1:4], time=t_exact[1:4])
  plt <- plt + geom_line(data = df_exact, aes(x=N,y=time), inherit.aes = F,
                         lty=2)  +
    geom_point(data = df_exact, aes(x=N,y=time), inherit.aes = F)
  return(plt)
}

