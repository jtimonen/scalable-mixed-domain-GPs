ker_eq <- function(x1, x2, z1, z2, ell) {
  ls <- 1 # sqrt(ell[z1] * ell[z2]) + 0.3
  dx <- (x1 / ell[z1] - x2 / ell[z2])^2 / ls
  exp(-0.5 * (dx))
}

# Non-separable kernel
nonsep_kernel <- function(x1, x2, z1, z2, ell = c(1, 0.7, 0.4)) {
  N1 <- length(x1)
  N2 <- length(x2)
  K <- matrix(0, N1, N2)
  for (i in 1:N1) {
    for (j in 1:N2) {
      K[i, j] <- ker_eq(x1[i], x2[j], z1[i], z2[j], ell)
    }
  }
  K
}

# True kernel function (nonseparable)
simulate_nonsep <- function(df) {
  mu0 <- rep(0, nrow(df))
  K <- nonsep_kernel(df$x, df$x, df$z, df$z)
  f <- MASS::mvrnorm(n = 1, mu0, K)
  fun <- data.frame(f)
  cbind(df, fun)
}

# Simulating data
simulate_input <- function() {
  x <- seq(-1, 1, by = 0.05)
  N <- length(x)
  x <- c(x, x, x) + 0.02 * rnorm(3 * N)
  z <- as.factor(rep(c(1, 2, 3), each = N))
  df <- data.frame(x = x, z = z, id = z)
  df <- as_tibble(df) %>% arrange(z, x)
}
simulate_obs <- function(df) {
  df$y <- df$f + rnorm(n = nrow(df), mean = 0, sd = 0.1)
  df
}
split_train_test <- function(df) {
  N <- nrow(df)
  n <- round(N / 2)
  idx <- sample(N, n)
  idx <- which(df$x > -0.3 & df$x < 0.5)
  list(
    train = df[setdiff(1:N, idx), ],
    test = df[idx, ]
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

plot_fit <- function(r, df, df_train, df_test) {
  c1 <- "gray30"
  r$plot(x_var = "x", color_by = NULL) +
    geom_point(
      data = df_train, mapping = aes(x = x, y = y),
      inherit.aes = FALSE, color = c1
    ) +
    geom_line(
      data = df, mapping = aes(x = x, y = y_pred_true), inherit.aes = F,
      color = "black"
    ) +
    geom_point(
      data = df_test, mapping = aes(x = x, y = y),
      inherit.aes = FALSE, pch = 4, color = c1
    ) + ylim(ymin, ymax) + theme(legend.position = "none")
}

gppred <- function(df, df_pred, sigma, delta) {
  K <- nonsep_kernel(df$x, df$x, df$z, df$z)
  Ks <- nonsep_kernel(df$x, df_pred$x, df$z, df_pred$z)
  Kss <- nonsep_kernel(df_pred$x, df_pred$x, df_pred$z, df_pred$z)
  K <- list(hi = K)
  Ks <- list(hi = t(Ks))
  Kss_diag <- list(hi = diag(Kss))
  yp <- lgpr:::fp_gaussian.compute(K, Ks, Kss_diag, sigma^2, delta, df$y)
  yp$mean[, 2]
}
