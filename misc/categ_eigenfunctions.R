library(rstan)

# Function
create_C_matrix_zs <- function(num_cat) {
  z <- seq_len(num_cat)
  lgpr:::STAN_kernel_const(z, z, 0, num_cat, get_stream())
}

# Eigendecomposition for zs kernel
decompose_zs <- function(n_cat) {
  evals <- rep(0.0, n_cat)
  evecs <- matrix(0.0, n_cat, n_cat)
  evals[1] <- 0.0
  evals[2:n_cat] <- n_cat/(n_cat-1.0)
  ev1 <- rep(1.0, n_cat)
  evecs[,1] <- ev1
  evecs[,2:n_cat] <- pracma::nullspace(t(ev1))
  list(
    values = evals,
    vectors = evecs
  )
}

# Check eigedecomposition correctness
check_decomp <- function(eig) {
  eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
}

NC <- 4
K <- create_C_matrix_zs(NC)
V <- eigen(K, symmetric = T)
V2 <- decompose_zs(NC)
print(check_decomp(V2))
