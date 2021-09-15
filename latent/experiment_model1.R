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
library(posterior)
library(cmdstanr)
rstan_options(javascript = FALSE)
rstan_options(auto_write = TRUE)

# Settings
N <- 100
# n_per_N <- 1000
# N <- 10
model_idx <- 1
chains <- 4
scale_bf <- 1.5
NUM_BF <- c(2, 4, 10, 30)
do_lgpr_marginal <- TRUE
backend <- "rstan"

# Simulate data using lgpr
# sd <- simulate_data(
#  N = N, t_data = seq(1, 5, length.out = n_per_N),
#  relevances = c(0, 1, 1),
#  covariates = c(2),
#  n_categs = c(2),
#  lengthscales = c(1.0, 1.0, 0.75), t_jitter = 0.2
# )
# dat <- sd@data

# Simulate
age <- seq(1, 5, length.out = N)
y <- sin(0.3 * age**2) + 0.2 * rnorm(N)
dat <- data.frame(age, y)
dat$y <- normalize_var(dat$y)

# Create model using lgpr
if (model_idx == 1) {
  form <- y ~ age
} else if (model_idx == 2) {
  form <- y ~ age + age | z
} else if (model_idx == 3) {
  form <- y ~ age + age | z + id
} else {
  form <- y ~ age + age | z + age | id
}
# prior <- list(ell = igam(4, 4))
prior <- list(ell = normal(0, 1))
model <- lgpr::create_model(form, dat, prior = prior, sample_f = TRUE)

# Exact fit(s)
N <- model@stan_input$num_obs
cat("N=", N, "\n", sep = "")
exact <- sample_exact(
  model,
  latent = FALSE, marginal = TRUE, backend = backend, refresh = 1000
)

# Approximate fits
NUM_BF <- c(3, 10, 30)
SCALE_BF <- c(1.5, 2.5)
K_plots <- list()
F_plots <- list()
PRES <- list()
conf_names <- c()
j <- 0
for (scale_bf in SCALE_BF) {
  j <- j + 1
  conf_names[j] <- paste0("scale_bf = ", scale_bf)

  # Sample approximate models
  approx <- sample_approx_alter_num_bf(model, NUM_BF, scale_bf,
    backend = backend, refresh = 1000
  )

  # Collect all fits
  fits <- c(approx$fits, exact)
  stan_dats <- approx$stan_dats

  # Results
  PRES[[j]] <- summarize_results(fits)
  K_plots[[j]] <- plot_kernelcomparison_eq(PRES[[j]]$p_means, stan_dats, 1)
  F_plots[[j]] <- plot_f_compare_separate(dat, fits, last_is_exact = TRUE)
}
names(K_PLOTS) <- conf_names
names(PRES) <- conf_names
names(FPLOTS) <- conf_names
