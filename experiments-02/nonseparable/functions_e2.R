ker_eq <- function(x1, x2, z1, z2, ell) {
  ls <- 1
  dx <- (x1 / ell[z1] - x2 / ell[z2])^2 / ls^2
  exp(-0.5 * (dx))
}

# Non-separable kernel
nonsep_kernel <- function(x1, x2, z1, z2, ell) {
  N1 <- length(x1)
  N2 <- length(x2)
  K <- matrix(0, N1, N2)
  for (i in 1:N1) {
    for (j in 1:N2) {
      K[i, j] <- ker_eq(x1[i], x2[j], z1[i], z2[j], ell)
    }
  }
  K * lgpr:::kernel_zerosum(z1, z2, 3)
}

# True kernel function (nonseparable)
simulate_nonsep <- function(df, ell = c(1, 0.7, 0.4)) {
  mu0 <- rep(0, nrow(df))
  K <- nonsep_kernel(df$x, df$x, df$z, df$z, ell)
  f <- MASS::mvrnorm(n = 1, mu0, K)
  fun <- data.frame(f)
  cbind(df, fun)
}

# Simulate linear data
simulate_fast <- function(df) {
  G <- length(unique(df$z))
  f <- (-0.5 + 0.2 * sin(6 * df$x * as.numeric(df$z) / G))
  f <- as.data.frame(f)
  cbind(df, f)
}

# Simulating data
simulate_input <- function() {
  x <- seq(-1, 1, by = 0.07)
  N <- length(x)
  x <- c(x, x, x) + 0.01 * rnorm(3 * N)
  z <- as.factor(rep(c(1, 2, 3), each = N))
  df <- data.frame(x = x, z = z, id = z)
  df <- as_tibble(df) %>% arrange(z, x)
}

# Simulating data
simulate_input_vary <- function(N, G) {
  x <- seq(-1, 1, length.out = N)
  x <- rep(x, G)
  z <- as.factor(rep(1:G, each = N))
  df <- data.frame(x = x, z = z, id = z)
  df <- as_tibble(df) %>% arrange(z, x)
}


simulate_obs <- function(df, sigma) {
  df$y <- df$f + rnorm(n = nrow(df), mean = 0, sd = sigma)
  df
}
split_train_test <- function(df, alt = TRUE, test_categ = c(1, 2, 3)) {
  N <- nrow(df)
  n <- round(N / 2)
  idx <- sample(N, n)
  if (alt) {
    idx <- which(df$x > -0.3 & df$x < 0.5 & df$z %in% test_categ)
  }
  is_test <- rep(FALSE, N)
  is_test[idx] <- TRUE
  df$is_test <- is_test
  list(
    train = df[setdiff(1:N, idx), ],
    test = df[idx, ],
    full = df
  )
}

create_dummy_x <- function(df, C = 0) {
  df$x1 <- df$x
  df$x2 <- df$x
  df$x3 <- df$x
  df$x1[which(df$z != 1)] <- C
  df$x2[which(df$z != 2)] <- C
  df$x3[which(df$z != 3)] <- C
  df
}

plot_fit <- function(r, df, df_train, df_test, has_true = TRUE, reverse = F) {
  c1 <- "firebrick"
  d <- r$quantiles_df()
  plt <- ggplot(d, aes(x = x, y = med, ymin = out_low, ymax = out_up)) +
    geom_ribbon(fill = "steelblue", col = "steelblue", alpha = 0.6) +
    geom_line(col = "steelblue") +
    geom_point(
      data = df_train, mapping = aes(x = x, y = y),
      inherit.aes = FALSE, color = c1
    )
  if (reverse) {
    plt <- ggplot(d, aes(x = x, y = med, ymin = out_low, ymax = out_up)) +
      geom_point(
        data = df_train, mapping = aes(x = x, y = y),
        inherit.aes = FALSE, color = c1
      ) +
      geom_ribbon(fill = "steelblue", col = "steelblue", alpha = 0.6) +
      geom_line(col = "steelblue")
  }
  if (has_true) {
    plt <- plt + geom_line(
      data = df, mapping = aes(x = x, y = y_pred_true), inherit.aes = F,
      color = "black"
    )
    plt <- plt + facet_grid(. ~ z)
  } else {
    plt <- plt + facet_wrap(. ~ z)
  }
  plt <- plt + geom_point(
    data = df_test, mapping = aes(x = x, y = y),
    inherit.aes = FALSE, pch = 4, color = c1
  ) +
    ylab("y") +
    theme(legend.position = "none")
  if (exists("ymin") & exists("ymax")) {
    plt <- plt + ylim(ymin, ymax)
  }
  plt
}

gppred <- function(df, df_pred, sigma, delta, ell) {
  K <- nonsep_kernel(df$x, df$x, df$z, df$z, ell)
  Ks <- nonsep_kernel(df$x, df_pred$x, df$z, df_pred$z, ell)
  Kss <- nonsep_kernel(df_pred$x, df_pred$x, df_pred$z, df_pred$z, ell)
  K <- list(hi = K)
  Ks <- list(hi = t(Ks))
  Kss_diag <- list(hi = diag(Kss))
  yp <- lgpr:::fp_gaussian.compute(K, Ks, Kss_diag, sigma^2, delta, df$y)
  yp$mean[, 2]
}

# Add fits to df
compute_rmse <- function(r) {
  r$quantiles_df() %>%
    mutate(sq_error = (med - y)^2) %>%
    group_by(is_test) %>%
    summarize(rmse = sqrt(mean(sq_error)))
}

compute_accuracy <- function(df, r1, r2) {
  e0 <- df %>%
    mutate(sq_error = (y_pred_true - y)^2) %>%
    group_by(is_test) %>%
    summarize(rmse = sqrt(mean(sq_error)))
  e0$kernel <- "TRUE"
  e1 <- compute_rmse(r1)
  e1$kernel <- "K1"
  e2 <- compute_rmse(r2)
  e2$kernel <- "K2"
  rbind(e0, e1, e2) %>% arrange(is_test)
}
