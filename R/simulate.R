# Generate data by drawing from exact GP
simulate_data <- function(N, N_indiv, sigma) {
  N_tp <- N / N_indiv
  ages <- c()
  for (idx in 1:N_indiv) {
    ages <- c(ages, sort(runif(N_tp, 0, 8)))
  }
  ids <- rep(1:N_indiv, each = N_tp)
  zs <- rep(1:3, each = N / 3)
  dat <- data.frame(as.factor(ids), ages, as.factor(zs))
  colnames(dat) <- c("id", "age", "z")
  mu0 <- rep(0.0, nrow(dat))
  K1 <- lgpr:::kernel_eq(dat$age, dat$age, ell = 2)
  K2_a <- lgpr:::kernel_zerosum(dat$z, dat$z, 3)
  K2_b <- lgpr:::kernel_eq(dat$age, dat$age, ell = 1)
  f1 <- MASS::mvrnorm(n = 1, mu0, K1)
  f2 <- MASS::mvrnorm(n = 1, mu0, K2_a * K2_b)
  effects <- list(f1 = f1, f2 = f2, f = f1 + f2)
  dat$y <- f1 + f2 + rnorm(length(f1), sd = sigma)
  list(
    data = dat,
    effects = effects
  )
}
