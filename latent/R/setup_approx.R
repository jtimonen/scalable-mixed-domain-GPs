# Create Stan input for bf model
setup_approx <- function(model, num_bf, scale_bf) {
  stopifnot(is(model, "lgpmodel"))
  si <- model@stan_input
  si_bf <- list(num_bf = num_bf, scale_bf = scale_bf)
  si <- c(si, si_bf)
  si_add <- stan_input_approx_precomp(si)
  c(si, si_add)
}

dollar <- function(x, field) {
  a <- x[[field]]
  if (is.null(a)) {
    stop(field, " was NULL")
  }
  a
}

# Create additional Stan input for approximate models
stan_input_approx_precomp <- function(stan_input) {
  num_bf <- dollar(stan_input, "num_bf")
  if (any(num_bf > 0)) {
    comps <- dollar(stan_input, "components")

    # Categorical decomposition stuff
    C_mats <- create_C_matrices(stan_input)
    C_decs <- decompose_C_matrices(C_mats)
    validate_stan_input_approx(C_decs, C_mats) # extra computation
    si_add <- list(
      C2_vecs = C_decs$vectors[[2]],
      C2_vals = C_decs$values[[2]],
      C3_vecs = C_decs$vectors[[3]],
      C3_vals = C_decs$values[[3]]
    )
    si_add_add <- list(
      C2_num = length(si_add$C2_vals),
      C3_num = length(si_add$C3_vals)
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
  VALS <- list()
  VECS <- list()
  SIZES <- rep(1, J)
  for (j in seq_len(J)) {
    eg <- eigen(C_matrices[[j]])
    VALS[[j]] <- dollar(eg, "values")
    VECS[[j]] <- dollar(eg, "vectors")
    SIZES[j] <- length(VALS[[j]])
  }
  list(values = VALS, vectors = VECS, sizes = SIZES)
}

# Validation
validate_stan_input_approx <- function(C_decs, C_mats) {
  C_vals <- dollar(C_decs, "values")
  C_vecs <- dollar(C_decs, "vectors")
  J <- length(C_mats)
  for (j in seq_len(J)) {
    cat("\n----------------- j =", j, "------------------\n")
    V <- C_vecs[[j]]
    D <- diag(C_vals[[j]])
    C_rec <- V %*% D %*% t(V)
    diff <- as.vector(C_mats[[j]] - C_rec)
    mae <- max(abs(diff))
    if (mae > 1e-12) {
      msg <- paste0("Numerical problem in decomposition ", j, ", MAE=", mae)
      msg <- paste0(msg, ". Please report a bug! (validate_stan_input_approx)")
      stop(msg)
    }
  }
}
