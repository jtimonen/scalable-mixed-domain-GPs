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
model <- create_model(y ~ age + age | z + age | id, dat)

# Create additional Stan input
M_bf <- 25
L_bf <- 3.0
stan_data <- setup_basisfun(model, M_bf = M_bf, L_bf = L_bf)
N <- stan_data$num_obs
cat("N =", N, "\n")

# Sample model
#m1 <- stan_model("stan/lgp_covariance.stan")
#f1 <- sampling(m1, stan_data, chains = CHAINS, refresh = 500)

# Predict exactly
st <- proc.time()
fp1 <-  fp_exact(model, f1)
t1 <- as.numeric((proc.time() - st)[3])
cat("Took \u001b[34;1m", t1, "\u001b[0m seconds.\n", sep="")

# Predict approximately
fp2 <- fp_approx(model, f1, stan_data)

