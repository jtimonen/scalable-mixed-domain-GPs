# Source all R files
for (f in dir("R")) {
  path <- file.path("R", f)
  source(path)
}

# Requirements
library(lgpr)
check_lgpr_version()
library(rstan)
library(ggplot2)
library(ggpubr)
library(posterior)
library(cmdstanr)
rstan_options(javascript = FALSE)
rstan_options(auto_write = TRUE)

# Settings
N <- 100
N_indiv <- 10
chains <- 4
NUM_BF <- c(4, 8, 16, 32, 64)
SCALE_BF <- c(1.2, 1.5, 2.5)
backend <- "cmdstanr" # "rstan"

# Simulate data using lgpr
sd <- simulate_data(
  N = N_indiv, t_data = seq(1, 5, length.out = N / N_indiv),
  relevances = c(0, 1, 1),
  covariates = c(2),
  n_categs = c(2),
  lengthscales = c(0.75, 1.0, 0.75), t_jitter = 0.2
)
dat <- sd@data
dat$y <- normalize_var(dat$y)

# Create model using lgpr
form <- y ~ age + age | z

# prior <- list(ell = igam(4, 4))
prior <- list(ell = normal(0, 1))
model <- lgpr::create_model(form, dat, prior = prior, sample_f = TRUE)

# Exact fit(s)
N <- model@stan_input$num_obs
cat("N=", N, "\n", sep = "")
#exact <- sample_exact(
#  model,
#  latent = FALSE, marginal = TRUE, backend = backend, refresh = 1000
#)

# Approximate fits
K_plots <- list()
F0_plots <- list()
F1_plots <- list()
F2_plots <- list()
aes0 <- aes(x = age, y = f_mean, group = z_x_fit, color = fit)
aes1 <- aes(x = age, y = f_mean, group = fit, color = fit)
aes2 <- aes(x = age, y = f_mean, group = z_x_fit, color = fit)

PRES <- list()
conf_names <- c()
j <- 0
for (scale_bf in SCALE_BF) {
  j <- j + 1
  conf_names[j] <- paste0("c = ", scale_bf)

  # Sample approximate models
  approx <- sample_approx_alter_num_bf(model, NUM_BF, scale_bf,
    backend = backend, refresh = 1000, adapt_delta=0.95)

  # Collect all fits
  fits <- c(approx$fits, exact)
  stan_dats <- approx$stan_dats

  # Results
  PRES[[j]] <- summarize_results(fits)
  K_plots[[j]] <- plot_kernelcomparison_eq(PRES[[j]]$p_means, stan_dats, 1)
  plt0 <- plot_f_compare_separate(dat, fits,
    last_is_exact = TRUE,
    aes = aes0, comp_idx = 0,
    plot_data = TRUE
  )
  plt1 <- plot_f_compare_separate(dat, fits,
    last_is_exact = TRUE,
    aes = aes1, comp_idx = 1
  )
  plt2 <- plot_f_compare_separate(dat, fits,
    last_is_exact = TRUE,
    aes = aes2, comp_idx = 2
  )
  F0_plots[[j]] <- lapply(plt0, function(x) {
    x + facet_wrap(. ~ z)
  })
  F1_plots[[j]] <- plt1
  F2_plots[[j]] <- lapply(plt2, function(x) {
    x + facet_wrap(. ~ z)
  })
}
names(K_plots) <- conf_names
names(F0_plots) <- conf_names
names(F1_plots) <- conf_names
names(F2_plots) <- conf_names
names(PRES) <- conf_names

# Final plots
pf1 <- ggarrange(plotlist = F0_plots[[1]], nrow = 1)
pf2 <- ggarrange(plotlist = F0_plots[[2]], nrow = 1)
pf3 <- ggarrange(plotlist = F0_plots[[3]], nrow = 1)
plt_f <- ggarrange(pf1, pf2, pf3,
  ncol = 1, labels = conf_names,
  label.x = -0.03, label.y = 1.00
)
plt_k <- ggarrange(plotlist = K_plots, nrow = 1, labels = conf_names)

plot_f_compare_separate(dat, fits, last_is_exact = TRUE, aes)[[1]]

