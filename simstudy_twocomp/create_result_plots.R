# Startup
source(normalizePath(file.path("..", "common.R")))
startup(create_dir = FALSE)
res_dir <- "res_to_210"

# Load results
N_TRAIN <- c(60, 90, 120, 150, 180, 210)
results <- list()
j <- 0
used_scales <- c(1.5, 2.5, 4.0)
used_nbfs <- c(4, 8, 16, 32, 64)
N_S <- length(used_scales)
N_B <- length(used_nbfs)
num_repl <- 30
num_chains <- 4
for (N_train in N_TRAIN) {
  j <- j + 1
  res_approx <- array(0.0, dim = c(4, num_repl, N_S, N_B))
  res_exact <- array(0.0, dim = c(4, num_repl))
  for (repl_idx in 1:num_repl) {
    ri_fn <- repl_idx
    # ri_fn <- 0
    fn <- file.path(res_dir, paste0("repl_", ri_fn, "/res_", N_train, ".rds"))
    res <- readRDS(fn)
    mpld <- res$rtables_test$mlpd
    total_time <- colSums(res$summary$runtimes)
    num_div <- res$summary$num_div
    rr <- rbind(total_time, num_div)
    formatted <- format_results(rr, used_scales, used_nbfs)

    res_approx[1, repl_idx, , ] <- 1 / num_chains * formatted$total_time$approx
    res_exact[1, repl_idx] <- 1 / num_chains * formatted$total_time$exact
    res_approx[2, repl_idx, , ] <- formatted$num_div$approx
    res_exact[2, repl_idx] <- formatted$num_div$exact
    res_approx[3, repl_idx, , ] <- res$rtables_train$mlpd$approx
    res_exact[3, repl_idx] <- res$rtables_train$mlpd$exact
    res_approx[4, repl_idx, , ] <- res$rtables_test$mlpd$approx
    res_exact[4, repl_idx] <- res$rtables_test$mlpd$exact
  }
  rownames(res_exact) <- c("time", "num_div", "mlpd_train", "mlpd_test")
  results[[j]] <- list(exact = res_exact, approx = res_approx)
}
names(results) <- N_TRAIN

# MLPD PLOTS
p_train <- plot_mlpd(results, used_scales, used_nbfs, train = TRUE)
p_test <- plot_mlpd(results, used_scales, used_nbfs, train = FALSE)
ggsave(p_train, file = "mlpd_train.pdf", width = 7.58, height = 5.05)
ggsave(p_test, file = "mlpd_test.pdf", width = 7.58, height = 5.05)

# RUNTIME PLOT
rt <- plot_runtimes(results, used_scales, used_nbfs, N_TRAIN, 1.5)
ggsave(rt, file = "runtimes-c1_5_to_210.pdf", width = 5.8, height = 2.5)
