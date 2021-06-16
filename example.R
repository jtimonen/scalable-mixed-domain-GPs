# Requirements
library(lgpr)
library(rstan)
rstan_options(javascript = FALSE)

# Source all R files
for(f in dir("R")) {
  path <- file.path("R", f)
  source(path)
}

# Create model
model <- create_model(y ~ age + sex, testdata_002)

# Create Stan input
stan_data <- setup_basisfun(model, M_bf = 20, c_bf = 1.5)
stan_model <- stan_model("stan/lgp_basisfun.stan")
fit <- sampling(stan_model, stan_data, iter =  1000, chains = 1)
