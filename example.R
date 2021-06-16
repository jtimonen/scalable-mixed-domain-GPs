# Requirements
library(lgpr)
library(rstan)

# Source all R files
for(f in dir("R")) {
  path <- file.path("R", f)
  source(path)
}

# Create model
model <- create_model(y ~ age + sex, testdata_002)

# Create Stan input
stan_data <- setup_basisfun(model, num_bf = 20, width = 3)
stan_model <- stan_model("lgp_covariance.stan")
fit <- sampling(stan_model, stan_data, iter =  1000, chains = 1)
