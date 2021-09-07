# Source all R files
for (f in dir("R")) {
  path <- file.path("R", f)
  source(path)
}

# Requirements
library(lgpr)
check_lgpr_version()
library(rstan)
library(ggplot2)
library(ggpubr)
rstan_options(javascript = FALSE)
rstan_options(auto_write = TRUE)


# Simulate data using lgpr
n_per_N <- 20
sd <- simulate_data(
  N = 10, t_data = seq(1, 5, length.out = n_per_N),
  relevances = c(0, 1, 1),
  covariates = c(2),
  n_categs = c(3),
  lengthscales = c(1.5, 0.75, 0.75), t_jitter = 0.2
)
dat <- sd@data

# Create model using lgpr
model <- create_model(y ~ age + age | z + id, dat, sample_f = TRUE)

# Create additional Stan input
num_bf <- 30
scale_bf <- 1.25
decs <- categorical_kernel_decompositions(model)
si_add <- additional_stan_input(model, num_bf, scale_bf, decs$decompositions)
stan_data <- c(model@stan_input, si_add)

# Test creating transformed data
expose_stanfuns()
phi_mats <- create_phi_mats(stan_data)
psi_mats <- create_psi_mats(stan_data, phi_mats)

# Create model and sample
sm1 <- stan_model("stan/lgp_latent_approx.stan")
f1 <- sampling(sm1, data = stan_data, cores = 4)

# Create model and sample
# sm2 <- stan_model("stan/lgp_latent_covariance.stan")
# f2 <- sampling(sm2, data = stan_data, cores = 4)

# Compare
N <- stan_data$num_obs
cat("N=", N, "\n", sep = "")
