# Startup
for (f in dir("R")) {
  source(file.path("R", f))
}
outdir <- startup("experiment3")
set.seed(3432) # for reproducibility of data simulation

# Settings
confs <- list()
j <- 0
for (num_bf in c(4, 8, 16, 32)) {
  j <- j + 1
  confs[[j]] <- create_configuration(num_bf, 1.5)
}

# Global setup
model_formula <- y ~ age + age | z
N <- 120
N_indiv <- 12 #
N_train <- N - N / N_indiv
N_test <- N - N_train
chains <- 4
refresh <- 1000

# Simulate data using lgpr
sd <- lgpr::simulate_data(
  N = N_indiv, t_data = seq(1, 5, length.out = N / N_indiv),
  relevances = c(0, 1, 1),
  covariates = c(2),
  n_categs = c(3),
  lengthscales = c(0.5, 0.3, 0.3), t_jitter = 0.3,
  snr = 10
)
dat <- sd@data
dat$y <- 100 + 10 * dat$y

# Split to train and test data
test_id <- 1
split <- lgpr:::split_by_factor(dat, test = test_id, var_name = "id")
train_dat <- split$train
test_dat <- split$test
N_train <- nrow(train_dat)
N_test <- nrow(test_dat)
cat("N_train=", N_train, ", N_test=", N_test, "\n", sep = "")
test_z <- test_dat$z[1]

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

# Compute test elpds and rmse
num_fits <- length(fits)
elpds_w1 <- rep(0.0, num_fits)
elpds_w2 <- rep(0.0, num_fits)
rmses <- rep(0.0, num_fits)
for (j in 1:num_fits) {
  elpds_w1[j] <- compute_elpd(fits[[j]], preds[[j]], y_star, 1)
  elpds_w2[j] <- compute_elpd(fits[[j]], preds[[j]], y_star, 2)
  rmses[j] <- compute_rmse(fits[[j]], preds[[j]], y_star)
}
res <- rbind(elpds_w1, elpds_w2, rmses)
colnames(res) <- names(fits)

# Plot denser predictions
arange <- range(dat$age)
x_dense <- lgpr::new_x(train_dat, seq(arange[1] - 0.3, arange[2] + 0.3, 0.1))
na_inds <- is.na(x_dense$id)
x_dense$id[na_inds] <- test_id
x_dense$z[na_inds] <- test_z
preds_dense <- compute_predictions(fits, x_dense)
plot_preds(train_dat, test_dat, preds_dense) + theme_bw() +
  theme(legend.position = "top")

# Runtimes plot
# rt <- plot_runtimes_wrt_N(PRES, NUM_BF, N_sizes, scale_bf)
# ggsave("res/exp3/times3.pdf", plot = rt, width = 5.5, height = 4.3)
