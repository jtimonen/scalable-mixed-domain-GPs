# Startup
backend <- "cmdstanr"
for (f in dir("R")) {
  source(file.path("R", f))
}
outdir <- startup("experiment1", backend)

# Settings
N <- 100
model_idx <- 1
chains <- 4
scale_bf <- 1.5
NUM_BF <- c(3, 5, 20, 50)
SCALE_BF <- c(1.1, 1.5, 2.5)
do_lgpr_marginal <- TRUE

# Simulate
age <- seq(1, 5, length.out = N)
y <- sin(0.3 * age**2) + 0.2 * rnorm(N)
dat <- data.frame(age, y)
dat$y <- normalize_var(dat$y)
dat$id <- as.factor(rep("ID1", length(y)))
# Create model using lgpr
form <- y ~ age
prior <- list(ell = normal(0, 1))
model <- lgpr::create_model(form, dat, prior = prior, sample_f = TRUE)

# Exact fit(s)
N <- model@stan_input$num_obs
cat("N=", N, "\n", sep = "")
exact <- sample_exact(
  model,
  latent = FALSE, marginal = TRUE, backend = backend, refresh = 1000
)

# Approximate fits
K_plots <- list()
F_plots <- list()
PRES <- list()
conf_names <- c()
j <- 0
for (scale_bf in SCALE_BF) {
  j <- j + 1
  conf_names[j] <- paste0("c = ", scale_bf)

  # Sample approximate models
  approx <- sample_approx_alter_num_bf(model, NUM_BF, scale_bf,
    backend = backend, refresh = 1000, adapt_delta = 0.95
  )

  # Collect all fits
  fits <- c(approx$fits, exact)
  stan_dats <- approx$stan_dats

  # Results
  PRES[[j]] <- summarize_results(fits)
  K_plots[[j]] <- plot_kernelcomparison_eq(PRES[[j]]$p_means, stan_dats, 1)
  F_plots[[j]] <- plot_f_compare_separate(dat, fits, last_is_exact = TRUE)
}
names(K_plots) <- conf_names
names(F_plots) <- conf_names
names(PRES) <- conf_names

# Final plots
pf1 <- ggarrange(plotlist = F_plots[[1]], nrow = 1)
pf2 <- ggarrange(plotlist = F_plots[[2]], nrow = 1)
pf3 <- ggarrange(plotlist = F_plots[[3]], nrow = 1)
plt_f <- ggarrange(pf1, pf2, pf3,
  ncol = 1, labels = conf_names,
  label.x = -0.03, label.y = 1.00
)
plt_k <- ggarrange(plotlist = K_plots, nrow = 1, labels = conf_names)
