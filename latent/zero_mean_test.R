# lgpr
library(lgpr)
if (packageVersion("lgpr") != "1.1.4") {
  stop("lgpr 1.1.4 is required!")
}

# Other requirements
library(rstan)
library(ggplot2)
library(ggpubr)
rstan_options(javascript = FALSE)
rstan_options(auto_write = TRUE)

source("R/functions.R")
expose_stanfuns()

N <- 100
B <- 10
PLOTS <- list()
SPLOTS <- list()
j <- 0
L <- 7
x <- seq(0, 5, length.out = N)
mat_B <- matrix(rep(seq_len(B), times = N), N, B, byrow = TRUE)
PHI <- STAN_basisfun_eq(x, mat_B, L)
LAM <- diag(exp(-1 + rnorm(B)))
PHI_zm <- PHI
for(b in 1:B) {
  PHI_zm[,b] <- PHI_zm[,b] - mean(PHI_zm[,b])
}
K <- PHI %*% LAM %*% t(PHI)
K_zm <- PHI_zm %*% LAM %*% t(PHI_zm)
