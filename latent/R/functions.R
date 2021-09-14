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


# UTILS -------------------------------------------------------------------

# Matrix rows to a list
matrix_to_list <- function(x) {
  m <- dim(x)[1]
  L <- list()
  for (i in seq_len(m)) {
    L[[i]] <- x[i, ]
  }
  return(L)
}


# STAN INPUT --------------------------------------------------------------

# Compute the C x C matrix for each term
create_C_matrices <- function(si) {
  STREAM <- get_stream()
  comp <- getElement(si, "components") # changes in 2.0
  Z <- getElement(si, "x_cat") # changes in 2.0
  Z_M <- getElement(si, "x_cat_num_levels") # changes in 2.0
  C_matrices <- list()
  J <- nrow(comp)
  iz <- comp[, 8] # changes in 2.0
  kz <- comp[, 2] # changes in 2.0
  for (j in seq_len(J)) {
    idx <- iz[j]
    if (idx != 0) {
      z <- unique(Z[idx, ])
      M <- Z_M[idx]
      C_j <- lgpr:::STAN_kernel_const(z, z, kz[j], M, STREAM)
    } else {
      C_j <- matrix(1.0, 1, 1)
    }
    C_matrices[[j]] <- C_j
  }
  return(C_matrices)
}

# Compute the eigendecompositions for the C x C matrices for each term
decompose_C_matrices <- function(C_matrices) {
  J <- length(C_matrices)
  cn <- names(C_matrices)
  VALS <- list()
  VECS <- list()
  SIZES <- rep(1, J)
  RANKS <- rep(1, J)
  for (j in seq_len(J)) {
    C_j <- C_matrices[[j]]
    if (!is.matrix(C_j)) {
      VALS[[j]] <- NA
      VECS[[j]] <- NA
      SIZES[j] <- 0
      RANKS[j] <- 0
    } else {
      eg <- eigen(C_matrices[[j]])
      evals <- getElement(eg, "values")
      inds <- which(abs(evals) > 1e-12) # remove "numerically zero" eigenvalues
      evals <- evals[inds]
      evecs <- getElement(eg, "vectors")[, inds, drop = FALSE]
      if (ncol(evecs) != length(evals)) {
        stop("error when creating categorical decomposition!")
      }
      VALS[[j]] <- evals
      VECS[[j]] <- evecs
      SIZES[j] <- nrow(evecs)
      RANKS[j] <- ncol(evecs)
    }
  }
  names(VALS) <- cn
  names(VECS) <- cn
  names(SIZES) <- cn
  names(RANKS) <- cn
  list(values = VALS, vectors = VECS, sizes = SIZES, ranks = RANKS)
}

# Validation
validate_C_decompositions <- function(C_decs, C_mats) {
  C_vals <- getElement(C_decs, "values")
  C_vecs <- getElement(C_decs, "vectors")
  C_ranks <- getElement(C_decs, "ranks")
  J <- length(C_mats)
  for (j in seq_len(J)) {
    if (C_ranks[j] > 0) {
      evals <- C_vals[[j]]
      V <- as.matrix(C_vecs[[j]])
      n <- length(evals)
      if (n == 1) {
        evals <- as.matrix(evals) # because diag(n) is not what we want
      }
      D <- diag(evals)
      C_rec <- V %*% D %*% t(V) # compute reconstruction
      diff <- as.vector(C_mats[[j]] - C_rec) # check if same as original
      mae <- max(abs(diff))
      if (mae > 1e-12) {
        msg <- paste0("Numerical problem in decomposition ", j, ", MAE=", mae)
        msg <- paste0(msg, ". Please report a bug! (validate_stan_input_approx)")
        stop(msg)
      }
    }
  }
}

# Compute categorical kernel decompositions
categorical_kernel_decompositions <- function(model) {
  check_model_compatibility(model)
  stan_input <- model@stan_input
  comps <- stan_input$components

  # Categorical decomposition stuff
  C_mats <- create_C_matrices(stan_input)
  names(C_mats) <- rownames(comps)
  C_decs <- decompose_C_matrices(C_mats)
  validate_C_decompositions(C_decs, C_mats) # extra computation, safety
  return(list(matrices = C_mats, decompositions = C_decs))
}

# Flatten arrays
unlist_C_decompositions <- function(C_decs) {
  C_vecs <- as.array(unlist(C_decs$vectors))
  C_vals <- as.array(unlist(C_decs$values))

  out <- list(
    C_vecs = as.array(C_vecs[which(!is.na(C_vecs))]),
    C_vals = as.array(C_vals[which(!is.na(C_vals))]),
    C_ranks = as.array(C_decs$ranks),
    C_sizes = as.array(C_decs$sizes)
  )
  out$C_rsp <- out$C_ranks * out$C_sizes
  return(out)
}

# Create additional Stan input needed for approx
additional_stan_input <- function(model, num_bf, scale_bf, decs) {
  check_model_compatibility(model)
  si_precomp <- unlist_C_decompositions(decs)
  comps <- model@stan_input$components
  type <- comps[, 1]
  halfrange <- function(x) {
    max(abs(x))
  }
  X_sd <- as.array(apply(model@stan_input$x_cont, 1, stats::sd))
  X_hr <- as.array(apply(model@stan_input$x_cont, 1, halfrange))
  num_xi <- pmax(as.numeric(type > 0) * num_bf, 1) * si_precomp$C_ranks
  si_bf <- list(
    num_bf = num_bf, scale_bf = scale_bf, X_sd = X_sd,
    X_hr = X_hr, num_xi = num_xi
  )
  si_add <- c(si_bf, si_precomp)
  return(si_add)
}


# FITTING STAN MODELS -----------------------------------------------------

# Sample approximate model
sample_approx <- function(model, num_bf, scale_bf, backend = "rstan", ...) {
  MODEL_FILE <- "stan/lgp_latent_approx.stan"
  decs <- categorical_kernel_decompositions(model)
  si_add <- additional_stan_input(model, num_bf, scale_bf, decs$decompositions)
  stan_data <- c(model@stan_input, si_add)
  if (backend == "cmdstanr") {
    sm <- cmdstanr::cmdstan_model(MODEL_FILE, include_paths = "stan")
    fit <- sm$sample(data = stan_data, ...)
  } else {
    sm <- rstan::stan_model(MODEL_FILE)
    fit <- rstan::sampling(sm, data = stan_data, ...)
  }
  list(
    stan_data = stan_data,
    fit = fit
  )
}

# Sample approximate model with various configurations
sample_approx_many <- function(model, num_bf, scale_bf,
                               backend = "rstan", ...) {
  stopifnot(is(model, "lgpmodel"))
  NUM_BF <- rep(num_bf, times = length(scale_bf))
  SCALE_BF <- rep(scale_bf, each = length(num_bf))
  J <- length(NUM_BF)
  fits <- list()
  stan_dats <- list()
  nams <- c()
  for (i in seq_len(J)) {
    num_bf <- NUM_BF[i]
    scale_bf <- SCALE_BF[i]
    conf_str <- paste0("* B = ", num_bf, ", c = ", scale_bf)
    cat(conf_str, "\n")
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
sample_exact <- function(model, latent = FALSE, marginal = TRUE, ...) {
  stopifnot(is(model, "lgpmodel"))
  fits <- list()
  nams <- c()
  if (latent) {
    sm <- rstan::stan_model("stan/lgp_latent.stan")
    fit <- rstan::sampling(sm, data = model@stan_input, ...)
    nams <- c(nams, "latent")
    fits <- c(fits, list(fit))
  }
  if (marginal) {
    fit <- lgpr::lgp(
      formula = model@model_formula,
      data = model@data, prior = model@full_prior
    )
    nams <- c(nams, "marginal")
    fits <- c(fits, list(fit))
  }
  names(fits) <- nams
  return(fits)
}


# WRAPPING STAN FUNCTIONS -------------------------------------------------

# Expose all Stan functions without creating a complete Stan model
expose_stanfuns <- function(STAN_HOME = "stan") {
  FILES <- c(
    "functions-utils.stan",
    "functions-kernels.stan",
    "functions-prior.stan",
    "functions-approx.stan"
  )

  # Create Stan model containing only a functions block with all the functions
  two_spaces <- "  "
  f_list <- lapply(file.path(STAN_HOME, FILES), FUN = readLines)
  functions <- paste(unlist(f_list), collapse = paste0("\n", two_spaces))
  functions <- paste0(two_spaces, functions)
  model_code <- paste(c("functions {", functions, "}"), collapse = "\n")
  header <- "// Automatically generated\n\n"
  model_code <- paste0(header, model_code, "\n")

  # Write Stan code to file and expose
  fn <- file.path(STAN_HOME, "all_functions.stan")
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

# Wraps STAN_build_f_latent
build_f_latent <- function(stan_data, PSI, alpha, ell, xi) {
  expose_stanfuns()
  N <- stan_data$num_obs
  seq_M <- STAN_seq_len(stan_data$num_bf)
  scale_bf <- stan_data$scale_bf
  X_hr <- as.array(stan_data$X_hr)
  num_xi <- as.array(stan_data$num_xi)
  comps <- matrix_to_list(stan_data$components)
  C_ranks <- as.array(stan_data$C_ranks)
  STAN_build_f_latent(
    N, comps, num_xi, C_ranks, seq_M, X_hr, scale_bf, PSI,
    alpha, ell, xi
  )
}

# Draw xi and build latent
draw_f_latent <- function(stan_data, PSI, alpha = NULL, ell = NULL) {
  if (is.null(alpha)) {
    alpha <- rep(1.0, stan_data$num_comps)
  }
  if (is.null(ell)) {
    ell <- rep(1.0, stan_data$num_ell)
  }
  alpha <- as.array(alpha)
  ell <- as.array(ell)
  P <- sum(stan_data$num_xi)
  xi <- rnorm(n = P)
  build_f_latent(stan_data, PSI, alpha, ell, xi)
}

# Approximate the EQ kernel
approximate_kernel_eq <- function(alpha, ell, stan_data, idx_x = 1) {
  td <- do_transformed_data(stan_data) # exposes stan functions
  PHI <- td$PHI_mats[[idx_x]]
  dj <- STAN_basisfun_eq_multipliers(alpha, ell, td$seq_B, td$L[idx_x])
  DELTA <- diag(dj)
  K <- PHI %*% DELTA %*% t(PHI)
  list(K_approx = K, x = stan_data$x_cont[idx_x, ])
}

# SUMMARIZING RESULTS --------------------------------------------------------

# Helper functions
get_times <- function(x) {
  if (is(x, "lgpfit")) x <- x@stan_fit
  if (is(x, "CmdStanMCMC")) {
    tims <- x$time()$chains$total
  } else {
    tims <- rowSums(get_elapsed_time(x))
  }
  return(tims)
}
t_mean <- function(x) stats::mean(get_times(x))
t_sd <- function(x) stats::sd(get_times(x))
get_ndiv <- function(x) sum(rstan::get_divergent_iterations(x))
get_pars <- function(x) {
  pars <- c("alpha", "ell", "sigma")
  if (is(x, "lgpfit")) x <- x@stan_fit
  if (is(x, "CmdStanMCMC")) {
    d <- x$draws(pars)
  } else {
    d <- posterior::as_draws(x, pars)
  }
  return(d)
}

# Get experiment results
summarize_results <- function(fits) {
  draws <- lapply(fits, get_pars)
  p_means <- sapply(draws, colMeans)
  colStds <- function(x) apply(x, 2, stats::sd)
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

# Compare approximate and exact EQ covariance functions
compare_kernels_eq <- function(pars_approx, pars, stan_data, idx_x = 1) {
  a <- approximate_kernel_eq(pars_approx[1], pars_approx[2], stan_data, idx_x)
  x <- stan_data$x_cont[idx_x, ]
  K <- STAN_kernel_eq(x, x, pars[1], pars[2])
  a$K <- K
  return(a)
}
