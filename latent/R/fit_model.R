
# Run sampling in either backend
run_sampling <- function(MODEL_FILE, stan_data, backend, ...) {
  if (backend == "cmdstanr") {
    sm <- cmdstanr::cmdstan_model(MODEL_FILE, include_paths = "stan")
    fit <- sm$sample(data = stan_data, ...)
  } else {
    sm <- rstan::stan_model(MODEL_FILE)
    fit <- rstan::sampling(sm, data = stan_data, ...)
  }
  return(fit)
}

# Sample approximate model
sample_approx <- function(exact_model, confs, backend, ...) {
  fits <- list()
  J <- length(confs)
  nams <- c()
  for (j in 1:J) {
    cf <- confs[[j]]
    stopifnot(isa(cf, "ExperimentConfiguration"))
    approx_model <- create_approx_model(exact_model, cf@num_bf, cf@scale_bf)
    si_add <- approx_model@add_stan_input
    stan_data <- c(exact_model@stan_input, si_add)
    fit <- run_sampling(approx_model@stan_file, stan_data, backend, ...)
    fits[[j]] <- new("ApproxModelFit",
      model = approx_model,
      fit = list(fit),
      backend = backend
    )
    nams <- c(nams, experiment_info(cf))
  }
  names(fits) <- nams
  return(fits)
}

# Create name for an approximate fit
create_fitname <- function(num_bf, scale_bf) {
  paste0("B = ", formatC(num_bf, width = 3))
}

# Sample approximate model with various configurations of num_bf
sample_approx_alter_num_bf <- function(model, NUM_BF, scale_bf,
                                       backend = "rstan", ...) {
  stopifnot(is(model, "lgpmodel"))
  J <- length(NUM_BF)
  fits <- list()
  stan_dats <- list()
  nams <- c()
  for (i in seq_len(J)) {
    num_bf <- NUM_BF[i]
    conf_str <- create_fitname(num_bf, scale_bf)
    cat("* ", conf_str, "\n", sep = "")
    sres <- sample_approx(model, num_bf, scale_bf, backend = backend, ...)
    fits[[i]] <- sres$fit
    stan_dats[[i]] <- sres$stan_data

    nams <- c(nams, conf_str)
  }
  names(fits) <- nams
  names(stan_dats) <- nams
  list(fits = fits, stan_dats = stan_dats)
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
    fit <- lgpr::lgp(
      formula = formula(model@model_formula@call),
      data = model@data, prior = model@full_prior, ...
    )
    nams <- c(nams, "marginal")
    fits <- c(fits, list(fit))
  }
  names(fits) <- nams
  return(fits)
}
