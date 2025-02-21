ker_eq <- function(x1, x2, z1, z2, ell) {
  ls <- 1 #sqrt(ell[z1] * ell[z2]) + 0.3
  dx <- (x1/ell[z1] - x2/ell[z2])^2 / ls
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
  x <- seq(-1, 1, by = 0.1)
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
create_dummy_x <- function(df, C = 0) {
  df$x1 <- df$x
  df$x2 <- df$x
  df$x3 <- df$x
  df$x1[which(df$z != 1)] <- C
  df$x2[which(df$z != 2)] <- C
  df$x3[which(df$z != 3)] <- C
  df
}
