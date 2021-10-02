# Startup
backend <- "cmdstanr"
for (f in dir("R")) {
  source(file.path("R", f))
}
outdir <- startup("experiment3", backend)

# Settings
confs <- list()
j <- 0
for (num_bf in c(8, 16, 32)) {
  j <- j + 1
  confs[[j]] <- create_configuration(num_bf, 1.5)
}

# Global setup
model_formula <- y ~ age + age | z
prior <- list(ell = normal(0, 1))
N_train <- 60
N_indiv <- N_train / 10
N_test <- 30
chains <- 4
N <- N_train + N_test

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
split <- lgpr:::split_random(dat, p_test = N_test / N)
train_dat <- split$train
test_dat <- split$test
N_train <- nrow(train_dat)
N_test <- nrow(test_dat)
cat("N_train=", N_train, ", N_test=", N_test, "\n", sep = "")

# Create model using lgpr
exact_model <- lgpr::create_model(model_formula, train_dat,
  prior = prior,
  sample_f = TRUE
)

# Fit exact model
if (N_train <= 200) {
  exact <- sample_exact(
    exact_model,
    latent = FALSE, marginal = TRUE, backend = backend, refresh = 1000
  )
}

results <- list()
j <- 0

# Fit approximate models
approx <- sample_approx(exact_model, confs,
  backend = backend, refresh = 1000,
  adapt_delta = 0.95
)

for (su in setups) {
  j <- j + 1

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

# Compute test elpd
num_fits <- length(fits)
ELPD <- rep(0.0, num_fits)
PRED <- list()
for (j in 1:num_fits) {
  fit <- fits[[j]]
  if (isa(fit, "lgpfit")) {
    num_bf <- NULL
  } else {
    num_bf <- NUM_BF[j]
  }
  elp <- compute_elpd(model, fit, test_dat, num_bf, scale_bf)
  ELPD[j] <- mean(elp$lpd)
  PRED[[j]] <- elp$pred
}
names(ELPD) <- names(fits)
names(PRED) <- names(fits)

x_test <- lgpr::new_x(train_dat, seq(0, 10, 0.2))

# Function for plotting
plot_pred_test <- function(train_dat, test_dat, PRED, cols) {
  N1 <- nrow(train_dat)
  N2 <- nrow(test_dat)
  dat <- rbind(train_dat, test_dat)
  dtype <- as.factor(c(rep("train", N1), rep("test", N2)))
  dat$dtype <- dtype
  plt <- ggplot2::ggplot(dat, aes(x = age, y = y, group = id, color = dtype))
  plt <- plt + geom_point() + facet_wrap(. ~ id)
  j <- 0
  for (p in PRED) {
    j <- j + 1
    if (isa(p, "Prediction")) {
      h <- colMeans(p@h)
      pdat <- cbind(p@x, h)
    } else {
      h <- colMeans(p@y_mean)
      pdat <- cbind(p@x, h)
    }
    plt <- plt + geom_line(
      data = pdat, aes(x = age, y = h, group = id),
      inherit.aes = FALSE, color = cols[j]
    )
  }
  return(plt)
}

plot_pred_test(
  train_dat, test_dat, PRED,
  c("orange", "blue", "purple", "black")
)
