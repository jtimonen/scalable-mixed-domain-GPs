
# Run sampling in either backend
run_sampling <- function(MODEL_FILE, stan_data, backend, ...) {
  if (backend == "cmdstanr") {
    stan_dir <- getOption("stan_dir")
    cat("model_file =", MODEL_FILE, "\n")
    cat("include_paths =", stan_dir, "\n")
    success <- FALSE
    tries <- 0
    while (success == FALSE && tries < 20) {
      tries <- tries + 1
      cat("* TRY", tries, ":\n")

      tryCatch(
        {
          sm <- cmdstanr::cmdstan_model(MODEL_FILE,
            include_paths = stan_dir,
            dir = stan_dir
          )
          fit <- sm$sample(data = stan_data, ...)
          success <- TRUE
        },
        error = function(e) {
          cat("* ERROR OCCURRED:\n")
          print(e)
          success <- FALSE
        }
      )
    }
  } else {
    tries <- 1
    sm <- rstan::stan_model(MODEL_FILE)
    fit <- rstan::sampling(sm, data = stan_data, ...)
  }
  return(list(fit = fit, tries = tries))
}

# Sample approximate model
sample_approx_beta <- function(exact_model, confs, dat, ...) {
  fits <- list()
  J <- length(confs)
  backend <- "cmdstanr"
  fits <- list()
  for (j in 1:J) {
    cf <- confs[[j]]
    print(cf)
    stopifnot(is(cf, "ExperimentConfiguration"))
    approx_model <- create_approx_beta_model(
      exact_model, cf@num_bf, cf@scale_bf, dat
    )
    stan_data <- get_full_stan_input(approx_model)
    sfit <- run_sampling(
      approx_model@stan_file, stan_data,
      backend, ...
    )
    fit <- sfit$fit

    fits[[j]] <- new("ApproxModelFit",
      model = approx_model,
      fit = list(fit),
      backend = backend,
      info = list(tries = sfit$tries)
    )
    names(fits)[j] <- experiment_info(cf)
  }
  return(fits)
}

# Create name for an approximate fit
create_fitname <- function(num_bf, scale_bf) {
  paste0("B = ", formatC(num_bf, width = 3))
}
