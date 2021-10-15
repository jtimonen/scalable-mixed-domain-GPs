#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  replication_idx <- 0
} else {
  replication_idx <- as.numeric(args[1])
}
cat("Currently in", getwd(), "\n")

# Startup
source(normalizePath(file.path("..", "common.R")))
outdir <- startup(replication_idx)

# Settings
confs <- list()
j <- 0
SCALES <- c(1.5, 2.5, 4)
NBFS <- c(4, 8, 16, 32, 64)
for (scale_bf in SCALES) {
  for (num_bf in NBFS) {
    j <- j + 1
    confs[[j]] <- create_configuration(num_bf, scale_bf)
  }
}

# Global setup
model_formula <- y ~ age + age | z
N_test <- 150
chains <- 4
iter <- 2000
refresh <- iter

N_train <- 120

# Generate data
sd <- simulate_data(N_train, N_test, 0.5)
train_dat <- sd$train_dat
test_dat <- sd$test_dat

# Create and fit exact model using lgpr
exact_model <- lgpr::create_model(model_formula, train_dat)
efit <- lgpr::sample_model(exact_model,
  chains = chains, refresh = refresh,
  iter = iter
)

# Fit approximate model with different configurations
afits <- sample_approx(exact_model, confs,
  refresh = refresh,
  chains = chains,
  adapt_delta = 0.95,
  iter_warmup = iter / 2,
  iter_sampling = iter / 2
)

# Summarize results
fits <- c(afits, list(efit))
names(fits)[length(fits)] <- "exact"
sumr <- summarize_results(fits)

# Predict
x_star <- create_x_dense(train_dat, test_dat)
preds_dense <- compute_predictions(fits, x_star)
y_star <- test_dat[["y"]]

# Plot
LLL <- length(preds_dense)
preds_to_plot <- preds_dense[c(1:4, LLL)]
plt <- plot_preds(train_dat, test_dat, preds_to_plot) + theme_bw() +
  theme(legend.position = "top") + ylab("y")
