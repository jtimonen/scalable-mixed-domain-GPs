# Create Stan input for bf model
setup_approx <- function(model, num_bf, scale_bf) {
  stopifnot(is(model, "lgpmodel"))
  si <- model@stan_input
  si_bf <- list(num_bf = num_bf, scale_bf = scale_bf)
  cn <- component_names(model)
  si_precomp <- stan_input_approx_precomp(c(si, si_bf), cn)
  c(si_bf, si_precomp)
}

dollar <- function(x, field) {
  a <- x[[field]]
  if (is.null(a)) {
    stop(field, " was NULL")
  }
  a
}

# Create additional Stan input for approximate models
stan_input_approx_precomp <- function(stan_input, comp_names) {
  num_bf <- dollar(stan_input, "num_bf")
  if (any(num_bf > 0)) {
    comps <- dollar(stan_input, "components")

    # Categorical decomposition stuff
    C_mats <- create_C_matrices(stan_input)
    names(C_mats) <- comp_names
    C_decs <- decompose_C_matrices(C_mats)
    validate_stan_input_approx(C_decs, C_mats) # extra computation
    si_add <- list(
      C2_vecs = t(C_decs$vectors[[2]]), # each row is one eigenvector
      C2_vals = as.array(C_decs$values[[2]]),
      C3_vecs = t(C_decs$vectors[[3]]), # each row is one eigenvector
      C3_vals = as.array(C_decs$values[[3]])
    )
    si_add_add <- list(
      C2_size = ncol(si_add$C2_vecs),
      C3_size = ncol(si_add$C3_vecs),
      C2_rank = nrow(si_add$C2_vecs),
      C3_rank = nrow(si_add$C3_vecs)
    )
    si_add <- c(si_add, si_add_add)
  } else {
    si_add <- c()
  }
  return(si_add)
}

# Compute the C x C matrix for each term
create_C_matrices <- function(si) {
  STREAM <- get_stream()
  comp <- dollar(si, "components") # changes in 2.0
  Z <- dollar(si, "x_cat") # changes in 2.0
  Z_M <- dollar(si, "x_cat_num_levels") # changes in 2.0
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
  for (j in seq_len(J)) {
    eg <- eigen(C_matrices[[j]])
    evals <- dollar(eg, "values")
    inds <- which(abs(evals) > 1e-12) # remove "numerically zero" eigenvalues
    evals <- evals[inds]
    evecs <- dollar(eg, "vectors")[, inds, drop = FALSE]
    if (ncol(evecs) != length(evals)) {
      stop("error when creating categorical decomposition!")
    }
    VALS[[j]] <- evals
    VECS[[j]] <- evecs
    SIZES[j] <- length(VALS[[j]])
  }
  names(VALS) <- cn
  names(VECS) <- cn
  names(SIZES) <- cn
  out <- list(values = VALS, vectors = VECS, sizes = SIZES)
  return(out)
}

# Validation
validate_stan_input_approx <- function(C_decs, C_mats) {
  C_vals <- dollar(C_decs, "values")
  C_vecs <- dollar(C_decs, "vectors")
  J <- length(C_mats)
  for (j in seq_len(J)) {
    V <- as.matrix(C_vecs[[j]])
    evals <- C_vals[[j]]
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
