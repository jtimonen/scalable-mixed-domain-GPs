# Startup
backend <- "cmdstanr"
for (f in dir("R")) {
  source(file.path("R", f))
}
outdir <- startup("experiment3", backend)

# Settings
N_sizes <- c(60) # , 80, 120, 160, 400) #, 1000, 2000)
chains <- 4
NUM_BF <- c(8, 16, 32)
scale_bf <- 1.5
j <- 0

PRES <- list()
conf_names <- c()

for (N in N_sizes) {
  conf_names[j] <- paste0("N = ", N)
  j <- j + 1
  N <- 2 * N
  N_indiv <- N / 10
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

  # Split to train and test data
  split <- lgpr:::split_random(dat, p_test = 10 / N)
  train_dat <- split$train
  test_dat <- split$test
  N_train <- nrow(train_dat)
  N_test <- nrow(test_dat)
  cat("N_train=", N_train, ", N_test=", N_test, "\n", sep = "")

  # Create model using lgpr
  form <- y ~ age + age | z

  # prior <- list(ell = igam(4, 4))
  prior <- list(ell = normal(0, 1))
  model <- lgpr::create_model(form, train_dat, prior = prior, sample_f = TRUE)

  # Exact fit(s)
  if (N_train <= 200) {
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
  fits <- c(approx$fits, exact)
  stan_dats <- approx$stan_dats

  # Results
  PRES[[j]] <- summarize_results(fits)
}

# names(PRES) <- conf_names

# Runtimes plot
# rt <- plot_runtimes_wrt_N(PRES, NUM_BF, N_sizes, scale_bf)
# ggsave("res/exp3/times3.pdf", plot = rt, width = 5.5, height = 4.3)


fit <- fits[[4]]
lpd <- compute_lpd(fit, df_star = test_dat)

fa <- fits[[1]]
fp <- pred_approx(model, fa, test_dat, NUM_BF[1], scale_bf)
# fp is a list (length S) of lists (length J) of vectors (length P)
