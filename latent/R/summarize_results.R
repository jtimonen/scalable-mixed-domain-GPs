# Helper functions
t_mean <- function(x) {
  mean(rowSums(get_elapsed_time(x)))
}
t_sd <- function(x) {
  stats::sd(rowSums(get_elapsed_time(x)))
}
get_pars <- function(x) {
  p <- extract(x, pars = c("alpha", "ell", "sigma"))
  return(cbind(p$alpha, p$ell, p$sigma))
}
get_ndiv <- function(x) {
  sum(rstan::get_divergent_iterations(x))
}
ensure_stanfit <- function(x) {
  if (isa(x, "lgpfit")) x <- x@stan_fit
  return(x)
}

# Get experiment results
summarize_results <- function(fits) {
  fits <- lapply(fits, ensure_stanfit)
  draws <- lapply(fits, get_pars)
  p_means <- sapply(draws, colMeans)
  colStds <- function(x) {
    apply(x, 2, stats::sd)
  }
  p_sds <- sapply(draws, colStds)
  list(
    draws = draws,
    p_means = p_means,
    p_sds = p_sds,
    t_means = sapply(fits, t_mean),
    t_sds = sapply(fits, t_sd),
    num_div = sapply(fits, get_ndiv)
  )
}
