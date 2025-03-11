library(lgpr2) # v0.0.3
library(lgpr)
library(tidyverse)
library(MASS)
library(ggpubr)
ggplot2::theme_set(ggplot2::theme_bw())
source("functions_e2.R")

sigma_true <- 0.1

set.seed(123)
m2 <- lgpr2:::LonModel$new(formula = y ~ gp(x, z) + gp(x))


run_exp <- function(N) {
  X <- simulate_input_vary(N, 10)
  df <- simulate_fast(X)
  df <- simulate_obs(df, sigma_true)
  df <- create_dummy_x(df)
  plt <- ggplot(df, aes(x = x, y = f, color = z)) +
    geom_line() +
    geom_point(mapping = aes(x = x, y = y, color = z)) +
    facet_grid(. ~ z)
  a <- split_train_test(df, test_categ = c(2, 3), alt = FALSE)
  df_train <- a$train
  df_test <- a$test
  df <- a$full

  # Fit
  t_start <- Sys.time()
  f2 <- m2$fit(
    data = df_train, chains = 1, iter_sampling = 4, iter_warmup = 4,
    refresh = 1
  )
  t_end <- Sys.time()
  t_dur <- as.numeric(t_end - t_start, units = "secs")
  t_stan <- f2$get_stan_fit()$time()$total

  # r2 <- f2$predict(df)$function_draws()$get_output()
  # r2 <- FunctionDraws$new(df, r2, "Kernel 1")

  # plot <- plot_fit(r2, df, df_train, df_test, FALSE, TRUE)
  size <- nrow(df_train)
  message(size)
  return(lst(t_dur, t_stan, size))
}
res <- NULL
sizes <- 2 * 10^c(1, 2, 3, 4, 5)
j <- 0
for (s in sizes) {
  j <- j + 1
  res[[j]] <- run_exp(s)
}
get_res <- function(x) {
  a <- c(x$size, x$t_dur, x$t_stan)
  df <- data.frame(a)
  df <- t(df)
  colnames(df) <- c("N", "t_full", "t_mcmc")
  t(df)
}
out <- data.frame(t(sapply(res, get_res)))
colnames(out) <- c("N", "t_full", "t_mcmc")
out$prop_overhead <- (out$t_full - out$t_mcmc) / (out$t_full)
out2 <- out[2:nrow(out), ]
plt_a <- out2 %>%
  ggplot(aes(x = N, y = t_mcmc)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  scale_x_log10() +
  ylab("MCMC time (seconds)")

plt_b <- out2 %>%
  ggplot(aes(x = N, y = prop_overhead)) +
  geom_line() +
  geom_point() +
  scale_x_log10() +
  ylab("Proportion of overhead time") +
  ylim(c(0.7, 1))

plt <- ggarrange(plt_a, plt_b, labels = c("a)", "b)"), nrow = 1)
ggsave(plt, file = "scaling_suppl.pdf", width = 8.1, height = 2.93)
