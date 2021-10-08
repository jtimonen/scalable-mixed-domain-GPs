# Startup
r_dir <- "../R/"
options(stan_dir = "../stan/")
for (f in dir(r_dir)) {
  source(file.path(r_dir, f))
}
source("plotting_01.R")
outdir <- startup()
set.seed(38832) # for reproducibility of data simulation

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
N <- 90
N_indiv <- 9 #
N_train <- N - N / N_indiv
N_test <- N - N_train
chains <- 4
iter <- 200
refresh <- iter

# Generate data
sd <- simulate_data(N, N_indiv, 0.3)
dat <- sd$data

# Split to train and test data
i_test <- which(dat$id %in% c(3, 6, 9))
split <- lgpr:::split_data(dat, i_test = i_test)
train_dat <- split$train
test_dat <- split$test
N_train <- nrow(train_dat)
N_test <- nrow(test_dat)
cat("N_train=", N_train, ", N_test=", N_test, "\n", sep = "")
test_z <- test_dat$z[1]

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
em <- compute_metrics(fits, preds, y_star, SCALES, NBFS)


# Plot denser predictions
arange <- range(dat$age)
xvals <- seq(arange[1] - 0.5, arange[2] + 0.5, 0.05)
age_dense <- rep(xvals, times = N_indiv)
id_dense <- rep(1:N_indiv, each = length(xvals))
z_dense <- rep(1:3, each = 3 * length(xvals))
x_dense <- data.frame(as.factor(id_dense), age_dense, as.factor(z_dense))
colnames(x_dense) <- c("id", "age", "z")
preds_dense <- compute_predictions(fits, x_dense)

# Plots
# pred_plots <- list()
# for (j in 1:length(apreds_dense)) {
#  pd <- c(apreds_dense[[j]], epred_dense)
#  names(pd) <- c(names(apreds_dense[[j]]), "exact")
#  plt_pred <- plot_preds(train_dat, test_dat, pd) + theme_bw() +
#    theme(legend.position = "top") + ylab("")
#  pred_plots[[j]] <- plt_pred
# }
# names(pred_plots) <- names(apreds_dense)


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
