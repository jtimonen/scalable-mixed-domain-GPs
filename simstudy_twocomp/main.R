#!/usr/bin/env Rscript

# Main R script for simstudy_twocomp
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

N_TRAIN <- c(60, 90, 120, 150, 180, 210)
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
  y_star <- test_dat[["y"]]

  # Compute metrics
  em_test <- compute_metrics(fits, preds_test, y_star)
  rtables_test <- format_results(em_test, SCALES, NBFS)
  em_train <- compute_metrics(fits, preds_train, train_dat$y)
  rtables_train <- format_results(em_train, SCALES, NBFS)
  print(rtables_train)
  print(rtables_test)

  # Save results to RDS file
  all_results <- list(
    train_dat = train_dat,
    test_dat = test_dat,
    summary = sumr,
    rtables_train = rtables_train,
    rtables_test = rtables_test,
    sess_info = sessionInfo()
  )
  fn <- file.path(outdir, paste0("res_", N_train, ".rds"))
  cat(" * Saving to", fn, "\n")
  saveRDS(all_results, file = fn)
}
