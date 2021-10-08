# Startup
r_dir <- "../R/"
options(stan_dir = "../stan/")
for (f in dir(r_dir)) {
  source(file.path(r_dir, f))
}
source("simulate_01.R")
source("plotting_01.R")
outdir <- startup()
set.seed(3922) # for reproducibility of data simulation

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
apreds <- compute_predictions(afits, test_dat)
epred <- pred_exact(efit, test_dat) # takes long
preds <- c(apreds, list(pe))
names(preds)[length(preds)] <- "exact"
y_star <- test_dat[["y"]]
em <- compute_metrics(fits, preds, y_star)
rtables <- format_results(em, SCALES, NBFS)

# Exact
# fe <- fits[["exact"]]
# pe <- preds[["exact"]]
# mu <- colMeans(pe@y_mean) # mean estimates
# y_std_e <- colMeans(pe@y_std) # variance estimates
# y_std_e <- sqrt(sig2)
# sig <- lgpr::get_draws(fe, pars="sigma")
# sig <- fe@model@var_scalings$y@scale * sig

# Approx
# fa <- fits[[8]]
# pa <- preds[[8]]
# sf <- get_cmdstanfit(fa)
# y_scl <- fa@model@exact_model@var_scalings$y
# s_draws <- as.vector(posterior::merge_chains(sf$draws("sigma")))
# s_draws <- s_draws * y_scl@scale
# h_mean <- colMeans(pa@h)
# h_std <- apply(pa@h, 2, stats::sd)
# y_std_a <- sqrt(h_std**2 + mean(s_draws)**2)
# lpd <- compute_lpd.sampled_gaussian(pa, y_star, s_draws)[["way2"]]
# plot(y_std_a, ylim=c(0,1.5))
# lines(y_std_e, col="red")


# Plot denser predictions
x_dense <- create_x_dense(train_dat, test_dat)
preds_dense <- compute_predictions(fits, x_dense)

# Plots
plt1 <- plot_against_exact(1:4, train_dat, test_dat, preds_dense)
plt2 <- plot_against_exact(5:8, train_dat, test_dat, preds_dense)
plt3 <- plot_against_exact(9:12, train_dat, test_dat, preds_dense)
plt4 <- plot_against_exact(13:16, train_dat, test_dat, preds_dense)

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
