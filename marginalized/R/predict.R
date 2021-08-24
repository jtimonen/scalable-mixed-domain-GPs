# Create kernel computer  from model and parameter fit
create_kc <- function(model, stan_fit, x = NULL, reduce = NULL) {
  m <- model
  sf <- stan_fit
  STREAM <- rstan::get_stream()
  kc <- lgpr:::create_kernel_computer(m, sf, x, reduce, NULL, FALSE, STREAM)
  return(kc)
}

# Exact function posterior
fp_exact <- function(model, stan_fit, x = NULL, reduce = NULL) {
  kc <- create_kc(model, stan_fit, x, reduce)
  sigma <- lgpr::get_draws(stan_fit, reduce = reduce, pars="sigma")
  s2 <- as.vector(sigma)^2
  y <- lgpr:::get_y(model, original = FALSE)
  fp <- lgpr:::fp_gaussian(kc, s2, y, verbose=TRUE)
  return(fp)
}


# Approximate function posterior
fp_approx <- function(model, stan_fit, stan_data, x = NULL, reduce = NULL) {
  
}
