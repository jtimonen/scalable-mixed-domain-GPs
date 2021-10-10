
# Run sampling in either backend
run_sampling <- function(MODEL_FILE, stan_data, backend, ...) {
  if (backend == "cmdstanr") {
    stan_dir <- getOption("stan_dir")
    sm <- cmdstanr::cmdstan_model(MODEL_FILE, include_paths = stan_dir)
    fit <- sm$sample(data = stan_data, ...)
  } else {
    sm <- rstan::stan_model(MODEL_FILE)
    fit <- rstan::sampling(sm, data = stan_data, ...)
  }
  return(fit)
}

# Sample approximate model
sample_approx <- function(exact_model, confs, ...) {
  fits <- list()
  J <- length(confs)
  backend <- "cmdstanr"
  fits <- list()
  for (j in 1:J) {
    cf <- confs[[j]]
    print(cf)
    stopifnot(is(cf, "ExperimentConfiguration"))
    approx_model <- create_approx_model(
      exact_model, cf@num_bf, cf@scale_bf
    )
    stan_data <- get_full_stan_input(approx_model)
    fit <- run_sampling(
      approx_model@stan_file, stan_data,
      backend, ...
    )
    fits[[j]] <- new("ApproxModelFit",
      model = approx_model,
      fit = list(fit),
      backend = backend
    )
    names(fits)[j] <- experiment_info(cf)
  }
  return(fits)
}

# Create name for an approximate fit
create_fitname <- function(num_bf, scale_bf) {
  paste0("B = ", formatC(num_bf, width = 3))
}

# Sample exact model(s) for comparison
sample_exact <- function(model, latent = FALSE, marginal = TRUE,
                         backend = "rstan", ...) {
  stopifnot(is(model, "lgpmodel"))
  fits <- list()
  nams <- c()
  if (latent) {
    fit <- run_sampling("stan/lgp_latent.stan", model@stan_input, backend, ...)
    nams <- c(nams, "latent")
    fits <- c(fits, list(fit))
  }
  if (marginal) {
    nams <- c(nams, "marginal")
    fits <- c(fits, list(fit))
  }
  names(fits) <- nams
  return(fits)
}
