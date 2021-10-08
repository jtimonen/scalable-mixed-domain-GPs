# It is easier to navigate this file in the Rstudio text editor,
# where you can see the names of different sections and functions and
# navigate to them easily. Click the selection on the lower left corner,
# which says (Top Level) when the cursor is here.

# REQUIREMENTS ------------------------------------------------------------

# Check if model has correct format
check_model_compatibility <- function(model) {
  stopifnot(is(model, "lgpmodel"))
  ver <- model@info$lgpr_version
  if (ver < "1.1.4" || ver >= "2.0.0") {
    stop(
      "model was created with lgpr ", ver, ", but should be created with ",
      "lgpr version number at least 1.1.4 and less than 2.0.0!"
    )
  }
  TRUE
}

# Check lgpr package version
check_lgpr_version <- function() {
  ver <- packageVersion("lgpr")
  if (ver < "1.1.4" || ver >= "2.0.0") {
    stop(
      "using lgpr version ", ver, ", but ",
      "lgpr version number at least 1.1.4 and less than 2.0.0 is needed!"
    )
  }
  TRUE
}

# Startup for all experiments
startup <- function() {
  library(lgpr)
  check_lgpr_version()
  library(ggplot2)
  library(ggpubr)
  library(posterior)
  library(rstan)
  library(RColorBrewer)
  rstan::rstan_options(javascript = FALSE)
  rstan::rstan_options(auto_write = TRUE)
  library(cmdstanr)
  outdir <- file.path("results")
  if (!dir.exists(outdir)) dir.create(outdir)
  return(outdir)
}

# MISC UTILS ----------------------------------------------------------------

# Matrix rows to a list
matrix_to_list <- function(x) {
  m <- dim(x)[1]
  L <- list()
  for (i in seq_len(m)) {
    L[[i]] <- as.numeric(x[i, ])
  }
  return(L)
}

# Matrix rows to a list of lists
matrix_to_list_of_lists <- function(x) {
  m <- dim(x)[1]
  L <- list()
  for (i in seq_len(m)) {
    L[[i]] <- as.list(x[i, ])
  }
  return(L)
}

# Standardize to zero mean and unit variance
normalize_var <- function(x) (x - mean(x)) / stats::sd(x)


# WRAPPING STAN FUNCTIONS -------------------------------------------------

# Expose all Stan functions without creating a complete Stan model
expose_stanfuns <- function() {
  stan_dir <- getOption("stan_dir")
  FILES <- c(
    file.path("chunks", "functions-utils.stan"),
    file.path("chunks", "functions-kernels.stan"),
    file.path("chunks", "functions-prior.stan"),
    file.path("chunks", "functions-approx.stan")
  )

  # Create Stan model containing only a functions block with all the functions
  two_spaces <- "  "
  f_list <- lapply(file.path(stan_dir, FILES), FUN = readLines)
  functions <- paste(unlist(f_list), collapse = paste0("\n", two_spaces))
  functions <- paste0(two_spaces, functions)
  model_code <- paste(c("functions {", functions, "}"), collapse = "\n")
  header <- "// Automatically generated\n\n"
  model_code <- paste0(header, model_code, "\n")

  # Write Stan code to file and expose
  fn <- file.path(stan_dir, "all_functions.stan")
  cat(model_code, file = fn)
  rstan::expose_stan_functions(fn, verbose = TRUE)
  file.remove(fn)
}

# Wraps STAN_create_psi_mats
create_psi_mats <- function(stan_data, PHI_mats) {
  N <- stan_data$num_obs
  M <- stan_data$num_bf
  num_xi <- as.array(stan_data$num_xi)
  comps <- matrix_to_list(stan_data$components)
  x_cat <- matrix_to_list(stan_data$x_cat)
  C_ranks <- as.array(stan_data$C_ranks)
  C_sizes <- as.array(stan_data$C_sizes)
  C_rsp <- as.array(stan_data$C_rsp)
  C_vals <- as.array(stan_data$C_vals)
  C_vecs <- as.array(stan_data$C_vecs)
  STAN_create_psi_mats(
    num_xi, PHI_mats, comps,
    x_cat, C_vals, C_vecs, C_ranks, C_sizes, C_rsp
  )
}

# Wraps computing all the transformed data
do_transformed_data <- function(stan_data) {
  expose_stanfuns()
  num_bf <- stan_data$num_bf
  num_obs <- stan_data$num_obs
  seq_B <- STAN_seq_len(num_bf)
  X_hr <- as.vector(stan_data$X_hr)
  mat_B <- matrix(rep(seq_B, each = num_obs), num_obs, num_bf, byrow = FALSE)
  L <- stan_data$scale_bf * X_hr
  x_cont <- matrix_to_list(stan_data$x_cont)
  PHI_mats <- STAN_create_basisfun_mats(x_cont, mat_B, L)
  PSI_mats <- create_psi_mats(stan_data, PHI_mats)
  list(
    seq_B = seq_B,
    mat_B = mat_B,
    L = L,
    PHI_mats = PHI_mats,
    PSI_mats = PSI_mats
  )
}

# SUMMARIZING RESULTS --------------------------------------------------------

# Helper functions
get_runtimes <- function(x) {
  if (is(x, "lgpfit")) x <- x@stan_fit
  if (is(x, "ApproxModelFit")) {
    fit <- get_cmdstanfit(x)
    tims <- fit$time()$chains$total
  } else {
    tims <- as.numeric(rowSums(get_elapsed_time(x)))
  }
  return(tims)
}
get_ndiv <- function(x) {
  if (is(x, "lgpfit")) x <- x@stan_fit
  if (is(x, "ApproxModelFit")) {
    fit <- x@fit[[1]]
    ndiv <- sum(fit$sampler_diagnostics()[, , "divergent__"])
  } else {
    ndiv <- sum(rstan::get_divergent_iterations(x))
  }
  return(ndiv)
}

# Get param draws
get_pars <- function(x) {
  pars <- c("alpha", "ell", "sigma")
  if (is(x, "lgpfit")) x <- x@stan_fit
  if (is(x, "ApproxModelFit")) {
    fit <- x@fit[[1]]
    d <- fit$draws(pars)
  } else {
    d <- posterior::subset_draws(posterior::as_draws(x), pars)
  }
  d <- posterior::merge_chains(d)
  return(d)
}

# Get experiment results
summarize_results <- function(fits) {
  draws <- lapply(fits, get_pars)
  p_means <- sapply(draws, function(x) apply(x, 3, mean))
  p_sds <- sapply(draws, function(x) apply(x, 3, stats::sd))
  colnames(p_means) <- names(fits)
  colnames(p_sds) <- names(fits)
  list(
    draws = draws,
    p_means = p_means,
    p_sds = p_sds,
    runtimes = sapply(fits, get_runtimes),
    num_div = sapply(fits, get_ndiv)
  )
}
