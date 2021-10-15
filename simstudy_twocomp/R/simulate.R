# Generate covariate data
simulate_data_x <- function(N_train, N_test) {
  N_train_tp <- N_train / 6
  N_test_tp <- N_test / 3
  df <- NULL
  for (idx in 1:9) {
    if (idx %in% c(7, 8, 9)) {
      k <- N_test_tp
      is_test <- rep(TRUE, k)
      ages <- sort(runif(k, 0, 10))
    } else {
      k <- N_train_tp
      is_test <- rep(FALSE, k)
      ages <- sort(runif(k, 0, 10))
    }
    ids <- rep(idx, k)
    z_val <- idx %% 3
    if (z_val == 0) {
      z_val <- 3
    }
    zs <- rep(z_val, k)
    df_idx <- data.frame(ids, zs, ages, is_test)
    df <- rbind(df, df_idx)
  }
  colnames(df) <- c("id", "z", "age", "is_test")
  df$id <- as.factor(df$id)
  df$z <- as.factor(df$z)
  df$is_test <- as.factor(df$is_test)
  return(df)
}

# Generate data by drawing from exact GP
simulate_data <- function(N_train, N_test, sigma) {
  dat <- simulate_data_x(N_train, N_test)
  mu0 <- rep(0.0, nrow(dat))
  K1 <- lgpr:::kernel_eq(dat$age, dat$age, ell = 2)
  K2_a <- lgpr:::kernel_zerosum(dat$z, dat$z, 3)
  K2_b <- lgpr:::kernel_eq(dat$age, dat$age, ell = 1)
  f1 <- MASS::mvrnorm(n = 1, mu0, K1)
  f2 <- MASS::mvrnorm(n = 1, mu0, K2_a * K2_b)
  effects <- list(f1 = f1, f2 = f2, f = f1 + f2)
  dat$y <- f1 + f2 + rnorm(length(f1), sd = sigma)
  dat$y <- 100 + 40 * dat$y

  # Split to train and test data
  i_test <- which(dat$is_test == TRUE)
  split <- lgpr:::split_data(dat, i_test = i_test)
  train_dat <- split$train
  test_dat <- split$test
  N_train <- nrow(train_dat)
  N_test <- nrow(test_dat)
  cat("N_train=", N_train, ", N_test=", N_test, "\n", sep = "")
  list(
    train_dat = train_dat,
    test_dat = test_dat,
    effects = effects
  )
}

# Denser set of prediction points for visualization
create_x_dense <- function(train_dat, test_dat, by = 0.2) {
  dat <- rbind(train_dat, test_dat)
  arange <- range(dat$age)
  xvals <- seq(arange[1] - 0.5, arange[2] + 0.5, by)
  age_dense <- rep(xvals, times = 9)
  id_dense <- rep(1:9, each = length(xvals))
  z_dense <- rep(c(1, 2, 3, 1, 2, 3, 1, 2, 3), each = length(xvals))
  x_dense <- data.frame(as.factor(id_dense), age_dense, as.factor(z_dense))
  colnames(x_dense) <- c("id", "age", "z")
  return(x_dense)
}
