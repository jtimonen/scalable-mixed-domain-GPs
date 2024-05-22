library(lgpr)
library(Matrix)

M <- 3
P <- 500
a <- simulate_data(N = M, t_data = seq(1, P, 1))
dat <- a@data
x <- dat$age
id <- dat$id
K1 <- lgpr:::kernel_eq(x, x, 3.2, 1.1)
K2 <- lgpr:::kernel_zerosum(id, id, M) * lgpr:::kernel_eq(x, x, 3.3, 1.2)
R1 <- rankMatrix(K1)
R2 <- rankMatrix(K2)
K <- K1 + K2
R <- rankMatrix(K)
