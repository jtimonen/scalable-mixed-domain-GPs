# Sample approximate model
sample_approx <- function(model, num_bf, scale_bf, ...) {
  decs <- categorical_kernel_decompositions(model)
  si_add <- additional_stan_input(model, num_bf, scale_bf, decs$decompositions)
  stan_data <- c(model@stan_input, si_add)
  # print(stan_data)
  sm <- stan_model("stan/lgp_latent_approx.stan")
  fit <- sampling(sm, data = stan_data, pars = "xi", include = FALSE, ...)
  list(
    stan_data = stan_data,
    fit = fit
  )
}
