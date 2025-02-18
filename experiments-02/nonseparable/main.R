library(lgpr2) # v0.0.3
library(lgpr)
library(tidyverse)
library(MASS)
library(ggpubr)
ggplot2::theme_set(ggplot2::theme_bw())

# Non-separable kernel
kfun <- function(x1, x2, z1, z2) {
  N1 <- length(x1)
  N2 <- length(x2)
  K <- matrix(0, N1, N2)
  for (i in 1:N1) {
    for (j in 1:N2) {
      if (z1[i] == z2[j]) {
        alpha <- 1
      } else {
        alpha <- 0.95
      }
      K[i, j] <- lgpr:::kernel_eq(x1[i], x2[j], alpha = alpha, ell = 0.5)
    }
  }
  K
}

# True kernel function (nonseparable)
nonsep_kernel2 <- function(df) {
  mu0 <- rep(0, nrow(df))
  K <- kfun(df$x, df$x, df$z, df$z)
  f <- MASS::mvrnorm(n = 1, mu0, K)
  fun <- data.frame(f)
  cbind(df, fun)
}

# True kernel function (nonseparable)
nonsep_kernel <- function(df, alpha_add) {
  K_separable <- lgpr:::kernel_eq(df$x, df$x, 1, 1) * lgpr:::kernel_zerosum(df$z, df$z, 3)
  K_add <- lgpr:::kernel_eq(df$x, df$x, alpha = alpha_add, ell = 0.5)
  mu0 <- rep(0, nrow(df))
  K <- K_separable + K_add
  f <- MASS::mvrnorm(n = 1, mu0, K)
  fun <- data.frame(f)
  cbind(df, fun)
}

# Simulating data
simulate_input <- function() {
  x <- seq(-1, 1, by = 0.1)
  N <- length(x)
  x <- c(x, x, x) + 0.02 * rnorm(3 * N)
  z <- as.factor(rep(c(1, 2, 3), each = N))
  df <- data.frame(x = x, z = z, id = z)
  df <- as_tibble(df) %>% arrange(z, x)
}

X <- simulate_input()
df <- nonsep_kernel2(X)
plt <- ggplot(df, aes(x = x, y = f, color = z)) +
  geom_line()
