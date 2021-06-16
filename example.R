# Requirements
library(lgpr)
library(cmdstanr)

# Source all R files
for (f in dir("R")) {
  path <- file.path("R", f)
  source(path)
}

# Create model
model <- create_model(y ~ age + age | sex, testdata_002)

# Create Stan input
stan_data <- setup_basisfun(model, M_bf = 20, L_bf = 3.0)

# Sample model 1
m1 <- cmdstan_model("stan/lgp_covariance.stan", include_paths = "stan")
f1 <- m1$sample(stan_data, chains = 1, refresh = 500)

# Sample model 2
m2 <- cmdstan_model("stan/lgp_basisfun.stan", include_paths = "stan")
f2 <- m2$sample(stan_data, chains = 1, refresh = 500)
