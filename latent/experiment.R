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

# Settings
n_per_N <- 10
N <- 8
model_idx <- 1
chains <- 4
scale_bf <- 3 / 2
NUM_BF <- c(20, 40, 60, 80)

# Simulate data using lgpr
sd <- simulate_data(
  N = N, t_data = seq(1, 5, length.out = n_per_N),
  relevances = c(0, 1, 1),
  covariates = c(2),
  n_categs = c(2),
  lengthscales = c(1.5, 1.0, 0.75), t_jitter = 0.2
)
dat <- sd@data
normalize_var <- function(x) (x - mean(x)) / stats::sd(x)
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
prior <- list(ell = igam(4, 4))
model <- create_model(form, dat, prior = prior, sample_f = TRUE)

# Approximate fits
NUM_CONF <- length(NUM_BF)
AFITS <- list()
for (i in seq_len(NUM_CONF)) {
  cat("\n================================================================\n")
  cat("i=", i, "\n", sep = "")
  res <- sample_approx(model, 30, scale_bf,
    chains = chains,
    refresh = 500
  )
  AFITS[[i]] <- res$fit
}

# Exact fit
N <- model@stan_input$num_obs
cat("N=", N, "\n", sep = "")
if (N <= 200) {
  sm_exact <- stan_model("stan/lgp_latent.stan")
  fit_exact <- sampling(sm_exact,
    data = model@stan_input, chains = chains,
    pars = "eta", include = FALSE, refresh = 500
  )
} else {
  fit_exact <- NULL
}

fit_lgpr <- lgp(formula = form, data = dat, prior = prior)

# Collect results
names(AFITS) <- NUM_BF
ex <- list(fit_exact, fit_lgpr@stan_fit)
names(ex) <- c("lgp_latent", "lgpr_marginal")
res <- list(approx_fits = AFITS, exact_fits = ex)

# Get experiment results
parse_results <- function(res) {
  t_mean <- function(x) {
    mean(rowSums(get_elapsed_time(x)))
  }
  t_sd <- function(x) {
    stats::sd(rowSums(get_elapsed_time(x)))
  }
  get_pars <- function(x) {
    p <- extract(x, pars = c("alpha", "ell", "sigma"))
    return(cbind(p$alpha, p$ell, p$sigma))
  }
  get_ndiv <- function(x) {
    sum(rstan::get_divergent_iterations(x))
  }
  nams <- c(names(res$approx_fits), names(res$exact_fits))
  ALL_FITS <- c(res$approx_fits, res$exact_fits)
  names(ALL_FITS) <- nams
  draws <- lapply(ALL_FITS, get_pars)
  p_means <- sapply(draws, colMeans)
  colStds <- function(x) {
    apply(x, 2, stats::sd)
  }
  p_sds <- sapply(draws, colStds)
  list(
    names = nams,
    t_means = sapply(ALL_FITS, t_mean),
    t_sds = sapply(ALL_FITS, t_sd),
    draws = draws,
    p_means = p_means,
    p_sds = p_sds,
    num_div = sapply(ALL_FITS, get_ndiv)
  )
}

pres <- parse_results(res)
