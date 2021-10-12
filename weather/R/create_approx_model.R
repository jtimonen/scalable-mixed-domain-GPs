# Compute the C x C matrix for each term
create_C_matrices <- function(si) {
  STREAM <- get_stream()
  comp <- getElement(si, "components")
  Z <- getElement(si, "x_cat")
  Z_M <- getElement(si, "x_cat_num_levels")
  C_matrices <- list()
  J <- nrow(comp)
  iz <- comp[, 8]
  kz <- comp[, 2]
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

# Create approximate model similar to exact model
create_approx_model <- function(model, num_bf, scale_bf, id_var) {
  stan_dir <- getOption("stan_dir")
  num_obs <- lgpr:::get_num_obs(model)
  om <- lgpr:::get_obs_model(model)
  stan_file <- file.path(stan_dir, "agp_approx_with_intercepts.stan")
  decs <- categorical_kernel_decompositions(model)
  si_add <- additional_stan_input(
    model, num_bf, scale_bf,
    decs$decompositions
  )

  # Add input related to intercepts
  z_names <- lgpr:::covariate_info(model)$categorical[, 1]
  id <- as.numeric(model@data[[id_var]])
  if (length(id) == 0) {
    stop("id_var not found in data!")
  }
  si <- lgpr:::get_stan_input(model)
  num_ids <- length(unique(id))
  idx_expand_intercepts <- rep(0, num_obs)
  for (uid in 1:num_ids) {
    inds <- which(id == uid)
    idx_expand_intercepts[inds] <- uid
  }
  si_add_add <- list(
    idx_expand_intercepts = idx_expand_intercepts,
    num_ids = num_ids
  )

  # Return
  new("ApproxModel",
    exact_model = model,
    added_stan_input = c(si_add, si_add_add),
    stan_file = stan_file,
    obs_model = om
  )
}
