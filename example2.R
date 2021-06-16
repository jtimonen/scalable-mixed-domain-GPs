# Requirements
library(lgpr)
library(rstan)
library(ggplot2)
library(ggpubr)

# Source all R files
for (f in dir("R")) {
  path <- file.path("R", f)
  source(path)
}

# Simulate data
sd <- simulate_data(
  N = 6, t_data = seq(1, 5, by = 0.05),
  relevances = c(0, 1, 1),
  covariates = c(2),
  lengthscales = c(1, 0.5, 0.5), t_jitter = 0.5
)
dat <- sd@data

# Create model
model <- create_model(y ~ age + age | z, dat)

# Create Stan input
M_bf <- 25
L_bf <- 3.0
stan_data <- setup_basisfun(model, M_bf = M_bf, L_bf = L_bf)
N <- stan_data$num_obs
cat("num_obs =", N, "\n")

# Sample model 1
#m1 <- stan_model("stan/lgp_covariance.stan")
#f1 <- sampling(m1, stan_data, chains = CHAINS)

# Sample model 2
m2 <- stan_model("stan/lgp_basisfun.stan")
f2 <- sampling(m2, stan_data, chains = CHAINS)

# Comparison plot
#ag_name <- paste0("M=", M_bf, ", L=", L_bf)
#plt <- plot_params_comparison(f1, f2, ag_name, N)

#ggsave("compare.pdf", plot=plt, width=9, height=5)
