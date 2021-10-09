# Startup
for (f in dir("R")) {
  source(file.path("R", f))
}
outdir <- startup("experiment3")
set.seed(93232) # for reproducibility of data simulation

# Settings
confs <- list()
j <- 0
for (num_bf in c(4, 8, 16, 32)) {
  j <- j + 1
  confs[[j]] <- create_configuration(num_bf, 1.5)
}

# Global setup
model_formula <- y ~ age + age | z
N_train <- 120
N_indiv <- 12 # N_train / 10
N_test <- 120
chains <- 4
refresh <- 1000
N <- N_train + N_test

# Simulate data using lgpr
sd <- lgpr::simulate_data(
  N = N_indiv, t_data = seq(1, 5, length.out = N / N_indiv),
  relevances = c(0, 1, 1),
  covariates = c(2),
  n_categs = c(3),
  lengthscales = c(0.5, 0.75, 0.75), t_jitter = 0.3,
  snr = 10
)
dat <- sd@data
dat$y <- 100 + 10 * dat$y

# Split to train and test data
split <- lgpr:::split_random(dat, p_test = N_test / N)
train_dat <- split$train
test_dat <- split$test
N_train <- nrow(train_dat)
N_test <- nrow(test_dat)
cat("N_train=", N_train, ", N_test=", N_test, "\n", sep = "")

# Create and fit exact model using lgpr
exact_model <- lgpr::create_model(model_formula, train_dat)
efit <- lgpr::sample_model(exact_model, chains = chains, refresh = refresh)

# Fit approximate model with different configurations
afits <- sample_approx(exact_model, confs,
  refresh = refresh,
  chains = chains,
  adapt_delta = 0.95
)

# Collect results
fits <- c(afits, list(efit))
names(fits)[length(fits)] <- "exact"
results <- summarize_results(fits)

# Predict
preds <- compute_predictions(fits, test_dat)
y_star <- test_dat[["y"]]

# Compute test elpds
num_fits <- length(fits)
elpds <- rep(0.0, num_fits)
for (j in 1:num_fits) {
  elpds[j] <- compute_lpd(fits[[j]], preds[[j]], y_star, 1)
}

print(elpds)

# Plot denser predictions
arange <- range(dat$age)
x_dense <- lgpr::new_x(train_dat, seq(arange[1] - 0.3, arange[2] + 0.3, 0.1))
preds_dense <- compute_predictions(fits, x_dense)
plot_preds(train_dat, test_dat, preds_dense)

# Runtimes plot
# rt <- plot_runtimes_wrt_N(PRES, NUM_BF, N_sizes, scale_bf)
# ggsave("res/exp3/times3.pdf", plot = rt, width = 5.5, height = 4.3)
