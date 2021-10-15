# Startup
source(normalizePath(file.path("..", "common.R")))
startup(create_dir = FALSE)

used_scales <- c(1.5, 2.5, 4.0)
used_nbfs <- c(4, 8, 16, 32, 64)
N_S <- length(used_scales)
N_B <- length(used_nbfs)
num_repl <- 30
num_chains <- 4

# Load results
res_dir <- "res_to_210"
N_TRAIN <- c(60, 90, 120, 150, 180, 210)
df <- NULL
for (N_train in N_TRAIN) {
  times <- array(0.0, dim = c(num_repl, N_B + 1))
  for (repl_idx in 1:num_repl) {
    ri_fn <- repl_idx
    fn <- file.path(res_dir, paste0("repl_", ri_fn, "/res_", N_train, ".rds"))
    res <- readRDS(fn)
    mpld <- res$rtables_test$mlpd
    total_time <- colSums(res$summary$runtimes)
    num_div <- res$summary$num_div
    rr <- rbind(total_time, num_div)
    formatted <- format_results(rr, used_scales, used_nbfs)
    ttt <- c(formatted$total_time$approx[1, ], formatted$total_time$exact)
    times[repl_idx, ] <- 1 / num_chains * ttt
  }
  model <- names(ttt)
  t_mean <- colMeans(times)
  t_std <- apply(times, 2, stats::sd)
  df_j <- data.frame(as.factor(model), t_mean, t_std)
  df_j$n_train <- rep(N_train, nrow(df_j))
  colnames(df_j)[1] <- "model"
  df <- rbind(df, df_j)
}


# Load more results
res_dir <- "res_to_4000"
N_TRAIN <- c(300, 500, 1000, 2000) # , 4000)
for (N_train in N_TRAIN) {
  times <- array(0.0, dim = c(num_repl, N_B))
  for (repl_idx in 1:num_repl) {
    ri_fn <- repl_idx
    fn <- file.path(res_dir, paste0("repl_", ri_fn, "/res_", N_train, ".rds"))
    res <- readRDS(fn)
    mpld <- res$rtables_test$mlpd
    total_time <- colSums(res$summary$runtimes)
    num_div <- res$summary$num_div
    rr <- rbind(total_time, num_div)
    formatted <- format_results(rr, used_scales, used_nbfs)
    ttt <- c(formatted$total_time$approx[1, ])
    times[repl_idx, ] <- 1 / num_chains * ttt
  }
  model <- names(ttt)
  t_mean <- colMeans(times)
  t_std <- apply(times, 2, stats::sd)
  df_j <- data.frame(as.factor(model), t_mean, t_std)
  df_j$n_train <- rep(N_train, nrow(df_j))
  colnames(df_j)[1] <- "model"
  df <- rbind(df, df_j)
}

# Plot
plt <- ggplot(df, aes(
  x = n_train, y = t_mean, ymin = t_mean - t_std,
  ymax = t_mean + t_std, group = model, color = model,
  pch = model
)) +
  geom_line() +
  geom_point()
