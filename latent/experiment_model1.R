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
rstan_options(javascript = FALSE)
rstan_options(auto_write = TRUE)

# Settings
N <- 200
# n_per_N <- 1000
# N <- 10
model_idx <- 1
chains <- 4
scale_bf <- 1.5
NUM_BF <- c(2, 4, 10, 30, 100, 200)
do_lgpr_marginal <- TRUE

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
# prior <- list(ell = igam(4, 4))
prior <- list(ell = normal(0, 1))
model <- create_model(form, dat, prior = prior, sample_f = TRUE)

# Approximate fits
J <- length(NUM_BF)
fits <- list()
stan_dats <- list()
for (i in seq_len(J)) {
  cat("\n================================================================\n")
  cat("i=", i, "\n", sep = "")
  sres <- sample_approx(model, NUM_BF[i], scale_bf,
    chains = chains, refresh = 100
  )
  fits[[i]] <- sres$fit
  stan_dats[[i]] <- sres$stan_data
}
names(fits) <- paste0("num_bf = ", NUM_BF)
names(stan_dats) <- paste0("num_bf = ", NUM_BF)

# Exact fit
N <- model@stan_input$num_obs
cat("N=", N, "\n", sep = "")
if (FALSE) {
  sm_exact <- stan_model("stan/lgp_latent.stan")
  fit_exact <- sampling(sm_exact,
    data = model@stan_input, chains = chains,
    pars = "eta", include = FALSE, refresh = 500
  )
} else {
  fit_exact <- NULL
}
if (do_lgpr_marginal) {
  fit_lgpr <- lgp(formula = form, data = dat, prior = prior)
  nam <- names(fits)
  fits <- c(fits, fit_lgpr)
  names(fits) <- c(nam, "lgpr_marginal")
}

# Results
pres <- summarize_results(fits)
rownames(pres$p_means) <- c("alpha", "ell", "sigma")
rownames(pres$p_sds) <- c("alpha", "ell", "sigma")

plt_same <- plot_f_compare_same(dat, fits)
plt_separate <- plot_f_compare_separate(dat, fits, last_is_exact = TRUE)

# Compare kernels
fit_name <- "num_bf = 2"
stan_data <- stan_dats[[fit_name]]
pars_lgpr <- as.vector(pres$p_means[, "lgpr_marginal"])
pars_approx <- as.vector(pres$p_means[, fit_name])
fit <- fits[[fit_name]]
expose_stanfuns()
plot_kernelcomparison_eq(pars_lgpr, pars_lgpr, stan_data, 1)
text(2, 1, paste0("exact (black) vs. ", fit_name, " (red) \n same params"))
