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
B <- 3
PLOTS <- list()
SPLOTS <- list()
j <- 0
LLL <- c(2, 5, 10)
for (j in 1:3) {
  L <- LLL[j]
  x <- seq(0, 5, length.out = N)
  mat_B <- matrix(rep(seq_len(B), times = N), N, B, byrow = TRUE)
  PHI <- STAN_basisfun_eq(x, mat_B, L)

  df <- data.frame(rep(x, B), as.vector(PHI), as.factor(as.vector(mat_B)))
  colnames(df) <- c("x", "phi", "b")
  PLOTS[[j]] <- ggplot(df, aes(x = x, y = phi, group = b, color = b)) +
    geom_line() +
    ylab("phi_b") + theme(legend.position = "top")
  phi_sum <- 1 * PHI[, 1] + 0.2 * PHI[, 2] + 0.6 * PHI[, 3]
  df2 <- data.frame(x, phi_sum)
  SPLOTS[[j]] <- ggplot(df2, aes(x = x, y = phi_sum)) +
    geom_line() +
    ylab("1.0 * phi_1 + 0.2 * phi_2 + 0.6 * phi_3")
}
labs <- paste("L =", LLL)
plt <- ggarrange(plotlist = c(PLOTS, SPLOTS), labels = c(labs, NA, NA, NA))
