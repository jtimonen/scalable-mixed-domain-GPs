# Requirements
library(lgpr)
library(rstan)
library(ggplot2)
library(ggpubr)
rstan_options(javascript = FALSE)
rstan_options(auto_write = TRUE)

# Source all R files
for (f in dir("R")) {
  path <- file.path("R", f)
  source(path)
}
CHAINS <- 4

# Create model using lgpr
model <- create_model(y ~ age + age | sex, testdata_002)

# Create additional Stan input
M_bf <- 25
L_bf <- 3.0
stan_data <- setup_basisfun(model, M_bf = M_bf, L_bf = L_bf)
N <- stan_data$num_obs

# Sample model 1
m1 <- stan_model("stan/lgp_covariance.stan")
f1 <- sampling(m1, stan_data, chains = CHAINS)

# Sample model 2
m2 <- stan_model("stan/lgp_basisfun.stan")
f2 <- sampling(m2, stan_data, chains = CHAINS)

# Comparison plot
ag_name <- paste0("M=", M_bf, ", L=", L_bf)
plot_params_comparison(f1, f2, ag_name, N)

