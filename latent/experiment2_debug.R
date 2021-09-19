# Startup
backend <- "cmdstanr"
for (f in dir("R")) {
  source(file.path("R", f))
}
outdir <- startup("experiment2", backend)

# Settings
N <- 120
N_indiv <- 10
chains <- 4
num_bf <- 32
scale_bf <- 1.5

# Simulate data using lgpr
sd <- simulate_data(
  N = N_indiv, t_data = seq(1, 5, length.out = N / N_indiv),
  relevances = c(0, 1, 1),
  covariates = c(2),
  n_categs = c(3),
  lengthscales = c(0.7, 1.0, 0.75), t_jitter = 0.2
)
dat <- sd@data
dat$y <- normalize_var(dat$y)

# Create model using lgpr
form <- y ~ age + age | z

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
aes0 <- aes(x = age, y = f_mean, group = z_x_fit, color = fit)
aes1 <- aes(x = age, y = f_mean, group = fit, color = fit)
aes2 <- aes(x = age, y = f_mean, group = z_x_fit, color = fit)


# Sample approximate models
approx <- sample_approx(model, num_bf, scale_bf, backend,
  refresh = 1000, adapt_delta = 0.95
)

# Collect all fits
fits <- list(approx$fit, exact)
names(fits) <- c("approx", "exact")
stan_dats <- approx$stan_data

# Results
PRES <- summarize_results(fits)
K_plot <- plot_kernelcomparison_eq(PRES$p_means, stan_dats, 1)
plt0 <- plot_f_compare_separate(dat, fits,
  last_is_exact = TRUE,
  aes = aes0, comp_idx = 0,
  plot_data = TRUE
)
plt1 <- plot_f_compare_separate(dat, fits,
  last_is_exact = TRUE,
  aes = aes1, comp_idx = 1
)
plt2 <- plot_f_compare_separate(dat, fits,
  last_is_exact = TRUE,
  aes = aes2, comp_idx = 2
)
F0_plots[[j]] <- lapply(plt0, function(x) {
  x + facet_wrap(. ~ z)
})
F1_plots[[j]] <- plt1
F2_plots[[j]] <- lapply(plt2, function(x) {
  x + facet_wrap(. ~ z)
})
