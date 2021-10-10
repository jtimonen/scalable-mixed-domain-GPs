# Startup
r_dir <- normalizePath("../R")
options(stan_dir = normalizePath("../stan"))
for (f in dir(r_dir)) {
  source(file.path(r_dir, f))
}
library(lgpr)
library(ggplot2)
library(ggpubr)

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
    fn <- paste0("res_01_triton/repl_", repl_idx, "/res_", N_train, ".rds")
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
  rownames(res_exact) <- c("time", "num_div", "mpld_train", "mpld_test")
  results[[j]] <- list(exact = res_exact, approx = res_approx)
}
names(results) <- N_TRAIN

# Plot results
plot_res <- function(res, var_idx, scales, bfs) {
  ra <- res$approx[var_idx,,,]
  ra_mean <- apply(ra, c(2,3), mean)
  ra_std <- apply(ra, c(2,3), stats::sd)
  re <- res$exact[var_idx,]
  vname <- rownames(res$exact)[var_idx]
  re_mean <- mean(re)
  re_std <- mean(re)
  df_val_mean <- as.vector(ra_mean)
  df_val_std <- as.vector(ra_std)
  df_c <- rep(scales, times=length(bfs))
  df_b <- rep(bfs, each=length(scales))
  df <- data.frame(as.numeric(df_b), as.factor(df_c), df_val_mean, df_val_std)
  colnames(df) <- c("B", "c", "valm", "vals")
  plt <- ggplot(df, aes(x=B, y=valm, group=c,color=c)) +  geom_line() +
    geom_point() + ylab(vname) + scale_x_continuous(breaks=bfs) +
    geom_hline(yintercept = re_mean,  lty=2) +
    theme(legend.position = "top") +
    ylim(-5.25, -4.55)
}

plots_test_mpld <- list()
for(j in 1:6) {
  plots_test_mpld[[j]] <- plot_res(results[[j]], 4, used_scales, used_nbfs)
}
all <- ggarrange(plotlist=plots_test_mpld, labels=paste("N_train =", N_TRAIN))
