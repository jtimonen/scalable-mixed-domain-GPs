# Requirements
library(lgpr)
library(cmdstanr)

# Source all R files
for(f in dir("R")) {
  path <- file.path("R", f)
  source(path)
}

# Create model
model <- create_model(y ~ age + sex, testdata_002)

# Create Stan input
stan_data <- setup_basisfun(model, M_bf = 20, c_bf = 1.5)
mod <- cmdstan_model("stan/lgp_basisfun.stan")
fit <- mod$sample(stan_data, chains = 1)
