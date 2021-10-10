#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  replication_idx <- 0
} else {
  replication_idx <- as.numeric(args[1])
}

# Startup
r_dir <- "../R/"
options(stan_dir = "../stan/")
for (f in dir(r_dir)) {
  source(file.path(r_dir, f))
}
source("simulate_01.R")
source("plotting_01.R")
source("predict_01.R")
outdir <- startup(replication_idx)
cat("\n * Results will be saved to:", outdir, "\n")

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

N_TRAIN <- c(60, 90) # , 120, 150, 180, 210)
for (N_train in N_TRAIN) {

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
  preds_train <- compute_predictions(fits, train_dat)
  preds_test <- compute_predictions(fits, test_dat)

  # Save results to RDS file
  all_results <- list(
    train_dat = train_dat,
    test_dat = test_dat,
    summary = sumr,
    preds_train = preds_train,
    preds_test = preds_test,
    fits = fits
  )
  fn <- file.path(outdir, paste0("res_", N_train, ".rds"))
  cat(" * Saving to", fn, "\n")
  saveRDS(all_results, file = fn)
}
