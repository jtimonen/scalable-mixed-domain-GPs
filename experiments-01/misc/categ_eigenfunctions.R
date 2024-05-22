library(rstan)

# Function
create_C_matrix_cs <- function(num_cat, rho) {
  (1 - rho) * diag(num_cat) + rho * matrix(1, num_cat, num_cat)
}

# Normalize matrix columns
normalize_columns <- function(m) {
  for (j in seq_len(ncol(m))) {
    m[, j] <- m[, j] / sqrt(sum(m[, j]**2))
  }
  return(m)
}
# Eigendecomposition for CS kernel (alpha=1)
decompose_cs <- function(n_cat, rho) {
  evals <- rep(0.0, n_cat)
  evecs <- matrix(0.0, n_cat, n_cat)
  evals[1] <- 1 + (n_cat - 1) * rho
  evals[2:n_cat] <- 1 - rho
  ev1 <- rep(1.0, n_cat)
  evecs[, 1] <- ev1
  H <- contr.helmert(n_cat)
  # H <- normalize_columns(H)
  evecs[, 2:n_cat] <- H # pracma::nullspace(t(rep(1.0, n_cat)))
  evecs <- normalize_columns(evecs)
  list(
    values = evals,
    vectors = evecs
  )
}

# Check eigedecomposition correctness
check_decomp <- function(eig) {
  eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
}

NC <- 5
rho <- rnorm(1, 0, 2)
K <- create_C_matrix_cs(NC, rho)
V1 <- eigen(K, symmetric = T)
V2 <- decompose_cs(NC, rho)
print(check_decomp(V1))
print(check_decomp(V2))
print(rho)
