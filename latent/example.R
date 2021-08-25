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


# Simulate data using lgpr
n_per_N <- 13
sd <- simulate_data(
  N = 6, t_data = seq(1, 5, length.out = n_per_N),
  relevances = c(0, 1, 1),
  covariates = c(2),
  lengthscales = c(1, 0.5, 0.5), t_jitter = 0.5
)
dat <- sd@data


# Create model using lgpr
model <- create_model(y ~ age + age | z + age | id, dat, sample_f = TRUE)

# Source all R files
for (f in dir("R")) {
  path <- file.path("R", f)
  source(path)
}
CHAINS <- 4

# Create additional Stan input
num_bf <- 25
scale_bf <- 1.5
stan_data <- setup_approx(model, num_bf = num_bf, scale_bf = scale_bf)

# Create model and sample
sm <- stan_model("stan/lgp_latent_basisfun.stan")
fit <- sampling(sm, data = stan_data, iter = 2)
