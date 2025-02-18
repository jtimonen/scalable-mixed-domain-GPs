library(lgpr2) # v0.0.3
library(lgpr)
library(tidyverse)
library(MASS)

# True kernel function (nonseparable)
true_kernel <- function(df) {
  K1 <- lgpr:::kernel_eq(df$x, df$x, 1, 1) * lgpr:::kernel_zerosum(df$z, df$z, 3)
  K2 <- lgpr:::kernel_eq(df$x, df$x, alpha = 0.5, ell = 0.5)
  mu0 <- rep(0, nrow(df))
  f1 <- MASS::mvrnorm(n = 1, mu0, K1)
  f2 <- MASS::mvrnorm(n = 1, mu0, K2)
  f <- f1 + f2
  fun <- data.frame(f1, f2, f)
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
df <- true_kernel(X)

plt <- ggplot(df, aes(x = x, y = f, color = z)) +
  geom_line()
