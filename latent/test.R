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
rstan_options(javascript = FALSE)
rstan_options(auto_write = TRUE)


# Simulate data using lgpr
n_per_N <- 20
sd <- simulate_data(
  N = 10, t_data = seq(1, 5, length.out = n_per_N),
  relevances = c(0, 1, 1),
  covariates = c(2),
  n_categs = c(2),
  lengthscales = c(1.5, 0.75, 0.75), t_jitter = 0.2
)
dat <- sd@data

# Create model using lgpr
model <- create_model(y ~ age + age | z + id, dat, sample_f = TRUE)

# Create additional Stan input
num_bf <- 30
scale_bf <- 1.25

res <- sample_approx(model, num_bf, scale_bf, chains = 2)
stan_data <- res$stan_data
f1 <- res$fit

# Create model and sample
# sm2 <- stan_model("stan/lgp_latent_covariance.stan")
# f2 <- sampling(sm2, data = stan_data, cores = 4)

N <- stan_data$num_obs
cat("N=", N, "\n", sep = "")

# Compare functions
f_draws1 <- extract(f1, pars = "f_latent")$f_latent
f_sum_draws <- apply(f_draws1, c(1, 3), sum)
fl_1 <- apply(f_draws1, c(2, 3), mean)
df1 <- data.frame(t(fl_1))
df1 <- cbind(df1, colMeans(f_sum_draws), apply(f_sum_draws, 2, stats::sd))
colnames(df1) <- c("f1", "f2", "f3", "f", "f_sd")
df1 <- cbind(dat, df1)

# f_draws2 <- extract(f2, pars = "f_latent")$f_latent
# fl_2 <- apply(f_draws2, c(2, 3), mean)
# df2 <- data.frame(t(fl_2))
# colnames(df2) <- c("f1")
# df2 <- cbind(dat, df2)

# type <- as.factor(rep(c("approx", "exact"), each = N))
# df <- rbind(df1, df2)
# df$type <- type

p1 <- ggplot(df1, aes(x = age, y = f1)) +
  geom_line()
p2 <- ggplot(df1, aes(x = age, y = f2, group = z, color = z)) +
  geom_line()
p3 <- ggplot(df1, aes(x = age, y = f3, group = id)) +
  geom_line() +
  facet_wrap(. ~ id)

p4 <- ggplot(df1, aes(x = age, y = f, group = id)) +
  geom_line() +
  geom_ribbon(
    mapping = aes(x = age, ymax = f + 2 * f_sd, ymin = f - 2 * f_sd, group = id),
    fill = "steelblue1", color = "steelblue1"
  ) +
  geom_point(mapping = aes(x = age, y = y, group = id)) +
  facet_wrap(. ~ id)
