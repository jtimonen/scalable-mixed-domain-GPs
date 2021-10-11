library(rstan)

# Function
create_C_matrix_zs <- function(num_cat) {
  z <- seq_len(num_cat)
  lgpr:::STAN_kernel_const(z, z, 0, num_cat, get_stream())
}

# Normalize matrix columns
normalize_columns <- function(m) {
  for (j in seq_len(ncol(m))) {
    m[, j] <- m[, j] / sqrt(sum(m[, j]**2))
  }
  return(m)
}
# Eigendecomposition for zs kernel
decompose_zs <- function(n_cat) {
  evals <- rep(0.0, n_cat)
  evecs <- matrix(0.0, n_cat, n_cat)
  evals[1] <- 0.0
  evals[2:n_cat] <- n_cat / (n_cat - 1.0)
  ev1 <- rep(1.0, n_cat)
  evecs[, 1] <- ev1
  H <- contr.helmert(n_cat)
  H <- normalize_columns(H)
  evecs[, 2:n_cat] <- H # pracma::nullspace(t(rep(1.0, n_cat)))
  list(
    values = evals,
    vectors = evecs
  )
}

# Check eigedecomposition correctness
check_decomp <- function(eig) {
  eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
}

NC <- 6
K <- create_C_matrix_zs(NC)
V <- eigen(K, symmetric = T)
V2 <- decompose_zs(NC)
print(check_decomp(V2))
