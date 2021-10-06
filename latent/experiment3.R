# Startup
for (f in dir("R")) {
  source(file.path("R", f))
}
outdir <- startup("experiment3")
set.seed(98342) # for reproducibility of data simulation

# More functions

# Create data frame for ggplot
create_pred_plot_df <- function(preds) {
  J <- length(preds)
  df <- NULL
  model <- c()
  for (j in 1:J) {
    p <- preds[[j]]
    if (isa(p, "Prediction")) {
      h <- colMeans(p@h)
    } else {
      h <- colMeans(p@y_mean)
    }
    df_add <- cbind(p@x, h)
    df <- rbind(df, df_add)
    model <- c(model, rep(names(preds)[j], nrow(df_add)))
  }
  df$model <- as.factor(model)
  df$id_x_model <- paste(df$id, df$model)
  return(df)
}

# Modify id factor for labeling
modify_id_label <- function(df, levels, labels) {
  df$id <- factor(df$id,
    levels = levels(df$id),
    labels = levels(labels)
  )
  return(df)
}

# Function for plotting predictions
plot_preds <- function(train_dat, test_dat, preds) {
  pdat <- create_pred_plot_df(preds)
  labels <- as.factor(paste(paste("id =", formatC(pdat$id, width = 2)),
    paste("z =", pdat$z),
    sep = ", "
  ))
  levels <- levels(labels)
  pdat <- modify_id_label(pdat, levels, labels)
  pdat_exact <- pdat[which(pdat$model == "exact"), ]
  pdat_approx <- pdat[which(pdat$model != "exact"), ]
  train_dat <- modify_id_label(train_dat, levels, labels)
  test_dat <- modify_id_label(test_dat, levels, labels)
  num_fits <- length(levels(pdat_approx$model))
  my_colors <- RColorBrewer::brewer.pal(num_fits, "PuBu")[2:num_fits]
  plt <- ggplot2::ggplot(pdat_approx, aes(
    x = age, y = h, group = model,
    color = model
  )) +
    geom_line(lwd = 0.8) +
    facet_wrap(. ~ id) +
    scale_color_manual(values = my_colors)
  plt <- plt + geom_line(
    data = pdat_exact, aes(x = age, y = h),
    color = "orange", lty = 1, inherit.aes = FALSE,
    lwd = 0.3
  )
  plt <- plt + geom_point(
    data = train_dat, aes(x = age, y = y), inherit.aes = FALSE,
    col = "black", pch = 20
  )
  plt <- plt + geom_point(
    data = test_dat, aes(x = age, y = y), inherit.aes = FALSE,
    col = "black", pch = 4
  )
  return(plt)
}

# Settings
confs <- list()
js <- 0
for (scale_bf in c(1.5, 2.5)) {
  js <- js + 1
  confs_s <- list()
  jb <- 0
  for (num_bf in c(6, 12, 24, 48)) {
    jb <- jb + 1
    confs_s[[jb]] <- create_configuration(num_bf, scale_bf)
  }
  confs[[js]] <- confs_s
}

# Global setup
model_formula <- y ~ age + age | z
N <- 120
N_indiv <- 12 #
N_train <- N - N / N_indiv
N_test <- N - N_train
chains <- 4
iter <- 100
refresh <- iter


# Simulate data using lgpr
sd <- lgpr::simulate_data(
  N = N_indiv, t_data = seq(1, 5, length.out = N / N_indiv),
  relevances = c(0, 1, 1),
  covariates = c(2),
  n_categs = c(3),
  lengthscales = c(0.6, 0.3, 0.3), t_jitter = 0.3,
  snr = 10
)
dat <- sd@data

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
efit <- lgpr::sample_model(exact_model,
  chains = chains, refresh = refresh,
  iter = iter
)

# Fit approximate model with different configurations
afits <- sample_approx(exact_model, confs,
  refresh = refresh,
  chains = chains,
  adapt_delta = 0.95,
  iter_warmup = iter,
  iter_sampling = iter
)

# Collect results
results <- lapply(afits, function(x) {
  get_results(x, efit)
})

# Predict
apreds <- lapply(afits, function(x) {
  compute_predictions(x, test_dat)
})
epred <- compute_predictions(list(efit), test_dat)
y_star <- test_dat[["y"]]

# Compute test elpds and rmse
num_fits <- length(afits)
em <- list()
for (j in 1:num_fits) {
  fs <- c(afits[[j]], list(efit))
  names(fs) <- c(names(afits[[j]]), "exact")
  ps <- c(apreds[[j]], epred)
  names(ps) <- c(names(apreds[[j]]), "exact")
  em[[j]] <- compute_metrics(fs, ps, y_star)
}
names(em) <- names(afits)


# Plot denser predictions
arange <- range(dat$age)
x_dense <- lgpr::new_x(train_dat, seq(arange[1] - 0.5, arange[2] + 0.5, 0.05))
na_inds <- is.na(x_dense$id)
x_dense$id[na_inds] <- test_id
x_dense$z[na_inds] <- test_z
apreds_dense <- lapply(afits, function(x) {
  compute_predictions(x, x_dense)
})
epred_dense <- compute_predictions(list(efit), x_dense)

# Plots
pred_plots <- list()
for (j in 1:length(apreds_dense)) {
  pd <- c(apreds_dense[[j]], epred_dense)
  names(pd) <- c(names(apreds_dense[[j]]), "exact")
  plt_pred <- plot_preds(train_dat, test_dat, pd) + theme_bw() +
    theme(legend.position = "top") + ylab("")
  pred_plots[[j]] <- plt_pred
}
names(pred_plots) <- names(apreds_dense)


# Save plot
j <- 0
for (pp in pred_plots) {
  j <- j + 1
  nam <- names(pred_plots)[j]
  np <- gsub("=", "_", gsub("[.]", "-", gsub(" ", "_", nam)))
  fn <- file.path(outdir, paste0("preds_", N, "_", np, ".pdf"))
  cat("Saving", fn, "\n")
  ggsave(fn, plot = pp, width = 6.7, height = 4.8)
}

# Runtimes plot
# rt <- plot_runtimes_wrt_N(PRES, NUM_BF, N_sizes, scale_bf)
# ggsave("res/exp3/times3.pdf", plot = rt, width = 5.5, height = 4.3)

# Save results in Rdata
res_to_save <- results[c("runtimes", "num_div")]
res_to_save$metrics <- em
