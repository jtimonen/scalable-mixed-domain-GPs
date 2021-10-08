# Startup
r_dir <- "../R/"
options(stan_dir = "../stan/")
for (f in dir(r_dir)) {
  source(file.path(r_dir, f))
}
source("simulate_01.R")
source("plotting_01.R")
outdir <- startup()
set.seed(28392) # for reproducibility of data simulation

# Settings
confs <- list()
j <- 0
SCALES <- c(1.5, 2.5, 5, 10)
NBFS <- c(4, 8, 16, 64)
for (scale_bf in SCALES) {
  for (num_bf in NBFS) {
    j <- j + 1
    confs[[j]] <- create_configuration(num_bf, scale_bf)
  }
}

# Global setup
model_formula <- y ~ age + age | z
N_train <- 90
N_test <- 150
chains <- 4
iter <- 200
refresh <- iter

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

# Collect results
fits <- c(afits, list(efit))
names(fits)[length(fits)] <- "exact"
results <- summarize_results(fits)

# Predict
preds <- compute_predictions(fits, test_dat)
y_star <- test_dat[["y"]]
em <- compute_metrics(fits, preds, y_star)
rtables <- format_results(em, SCALES, NBFS)


# Plot denser predictions
x_dense <- create_x_dense(train_dat, test_dat)
preds_dense <- compute_predictions(fits, x_dense)

# Plots
pd <- preds_dense[c(1:4, length(preds_dense))]
plt_pred <- plot_preds(train_dat, test_dat, pd) + theme_bw() +
  theme(legend.position = "top") + ylab("")


# Save plot
j <- 0
for (pp in pred_plots) {
  j <- j + 1
  nam <- names(pred_plots)[j]
  np <- gsub("=", "_", gsub("[.]", "-", gsub(" ", "_", nam)))
  fn <- file.path(outdir, paste0("preds_", N, "_", np, ".pdf"))
  cat("Saving", fn, "\n")
  ggsave(fn, plot = pp, width = 5.28, height = 4.75)
}

# Runtimes plot
# rt <- plot_runtimes_wrt_N(PRES, NUM_BF, N_sizes, scale_bf)
# ggsave("res/exp3/times3.pdf", plot = rt, width = 5.5, height = 4.3)

# Save results in Rdata
res_to_save <- results[c("runtimes", "num_div")]
res_to_save$metrics <- em

# Create latex table
library(xtable)
xtable(em[[1]], digits = 4)
xtable(em[[2]], digits = 4)
