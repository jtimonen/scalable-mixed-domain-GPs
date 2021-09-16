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
N_indiv_sizes <- c(8, 12, 16, 20, 24)
N_sizes <- rep(240, length(N_indiv_sizes))
chains <- 4
NUM_BF <- c(8, 16, 32)
scale_bf <- 1.5
backend <- "cmdstanr" # "rstan"
j <- 0

PRES <- list()
conf_names <- c()

for (N in N_sizes) {
  conf_names[j] <- paste0("N = ", N)
  j <- j + 1
  N_indiv <- N_indiv_sizes[j]
  # Simulate data using lgpr
  sd <- simulate_data(
    N = N_indiv, t_data = seq(1, 5, length.out = N / N_indiv),
    relevances = c(0, 1, 1),
    covariates = c(2),
    n_categs = c(3),
    lengthscales = c(0.5, 0.75, 0.75), t_jitter = 0.2
  )
  dat <- sd@data
  dat$y <- normalize_var(dat$y)

  # Create model using lgpr
  form <- y ~ age + age | z + age | id

  # prior <- list(ell = igam(4, 4))
  prior <- list(ell = normal(0, 1))
  model <- lgpr::create_model(form, dat, prior = prior, sample_f = TRUE)

  # Exact fit(s)
  N <- model@stan_input$num_obs
  cat("N=", N, "\n", sep = "")
  if (N <= 0) {
    exact <- sample_exact(
      model,
      latent = FALSE, marginal = TRUE, backend = backend, refresh = 1000
    )
  }

  # Sample approximate models
  approx <- sample_approx_alter_num_bf(model, NUM_BF, scale_bf,
    backend = backend, refresh = 1000, adapt_delta = 0.95
  )

  # Collect all fits
  fits <- approx$fits # c(approx$fits, exact)
  stan_dats <- approx$stan_dats

  # Results
  PRES[[j]] <- summarize_results(fits)
}

names(PRES) <- paste0("N_indiv = ", N_indiv_sizes)

# Runtimes plot
rt <- plot_runtimes_wrt_N_indiv(PRES, NUM_BF, N_indiv_sizes, N, scale_bf)
ggsave("res/exp4/times4.pdf", plot = rt, width = 5.5, height = 4.3)

# Divergences
divs <- sapply(PRES, function(x) {
  getElement(x, "num_div")
})
